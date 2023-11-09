from sim.healthcare.diagnosis import *
from sim.healthcare.system import *
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tswab', 'get_intv_poc', 'get_perfect_dx', 'get_perfect_txi']


def fillup(fr, to, tar):
    if fr + to < tar:
        return 0, fr + to
    elif to < tar:
        return fr + to - tar, tar
    else:
        return fr, to


def get_intv_tswab(p, p_loss_swab=0.05):
    sys0 = p['system']
    alg1, alg2, alg3 = sys0.Public.Algorithms
    _, alg4 = sys0.Engaged.Algorithms

    swab = Specimen('Swab', p_loss_swab)
    ssm = alg1.SSM
    xpert = Test('Xpert', alg3.Xpert.Sens, alg3.Xpert.Spec, swab)
    xpert_sn = Test('Xpert_ss-', alg1.Xpert.Sens, alg1.Xpert.Spec, swab)

    cdx_bn = alg1.CDx

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_sn, cdx=cdx_bn)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx_bn)

    ent_pub = sys0.Public.Entry.copy()
    ent_eng = sys0.Engaged.Entry.copy()

    ent = sys0.Entry.copy()

    public = Sector(ent_pub, [alg1, alg2, alg3])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    p_txi = p['p_txi']

    return {
        'sys': sys,
        'p_txi': p_txi
    }


def get_intv_poc(p, target=0.8, p_txi_poc=0.95, p_loss_swab=0.05):
    sys0 = p['system']
    alg1, alg2, alg3 = sys0.Public.Algorithms
    _, alg4 = sys0.Engaged.Algorithms

    swab = Specimen('Swab', p_loss_swab)
    ssm = alg1.SSM
    xpert = Test('Xpert', alg3.Xpert.Sens, alg3.Xpert.Spec, swab)
    xpert_sn = Test('Xpert_ss-', alg1.Xpert.Sens, alg1.Xpert.Spec, swab)

    cdx_bn = alg1.CDx

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_sn, cdx=cdx_bn)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx_bn)

    ent_pub = sys0.Public.Entry.copy()

    p_sx, p_s, p_x = ent_pub

    p_s, p_x = fillup(p_s, p_x, target)
    p_sx, p_x = fillup(p_sx, p_x, target)
    ent_pub = np.array([p_sx, p_s, p_x])

    ent_eng = p['system'].Engaged.Entry.copy()
    ent_eng[1], ent_eng[0] = fillup(ent_eng[1], ent_eng[0], target)

    ent = p['system'].Entry.copy()

    public = Sector(ent_pub, [alg1, alg2, alg3])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    p_txi_pub, p_txi_eng, p_txi_pri = p['p_txi']

    if p_txi_poc is not None:
        dx = np.array([alg.dx(1, 0).TruePos for alg in public.Algorithms])
        dx *= ent_pub
        dx /= dx.sum()
        p_txi_pub = (dx * np.array([p_txi_poc, p_txi_pub, p_txi_poc])).sum()

        dx = np.array([alg.dx(1, 0).TruePos for alg in engaged.Algorithms])
        dx *= ent_eng
        dx /= dx.sum()
        p_txi_eng = (dx * np.array([p_txi_poc, p_txi_eng])).sum()

    txi = np.array([p_txi_pub, p_txi_eng, p_txi_pri])

    return {
        'sys': sys,
        'p_txi': txi
    }


def get_perfect_dx(p):
    cdx = Test('CDx', 1, 1)

    alg4 = Algorithm('CDx', cdx=cdx)

    sys0 = p['system']
    ent = sys0.Entry.copy()

    public = Sector(np.array([1]), [alg4])
    engaged = Sector(np.array([1]), [alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    return {
        'sys': sys,
        'p_txi': p['p_txi']
    }


def get_perfect_txi(p):
    return {
        'sys': p['system'],
        'p_txi': np.ones(3)
    }
