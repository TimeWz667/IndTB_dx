from sim.healthcare.diagnosis import *
from sim.healthcare.system import *
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tswab', 'get_intv_poc', 'get_perfect_dx', 'get_perfect_txi']


def fillup(fr, to, tar):
    if fr + to < tar:
        return 0, fr + to
    elif to < tar:
        return fr + to - tar , tar
    else:
        return fr, to


def get_intv_tswab(p):
    sputum = Specimen('Sputum', p['p_loss_sputum'])
    swab = Specimen('Swab', p['p_loss_swab'])
    ssm = Test('SSM', p['sens_ssm'], p['spec_ssm'], sputum)
    xpert = Test('Xpert', p['sens_xpert'], p['spec_xpert'], swab)
    xpert_ssm = Test('Xpert_ss-', p['sens_xpert_ss-'], p['spec_xpert'], swab)

    cdx = Test('CDx', p['sens_cdx'], p['spec_cdx'])

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_ssm, cdx=cdx)
    alg2 = Algorithm('SSM > CDx', ssm=ssm, cdx=cdx)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx)
    alg4 = Algorithm('CDx', cdx=cdx)

    ent_pub = np.array([
        p['p_ava_ssm_pub'] * p['p_ava_xpert_pub'],
        p['p_ava_ssm_pub'] * (1 - p['p_ava_xpert_pub']),
        (1 - p['p_ava_ssm_pub']) * p['p_ava_xpert_pub'],
        (1 - p['p_ava_ssm_pub']) * (1 - p['p_ava_xpert_pub'])
    ])

    ent_eng = np.array([
        p['p_ava_xpert_eng'],
        (1 - p['p_ava_xpert_eng'])
    ])

    ent = np.array([
        p['p_csi_pub'],
        (1 - p['p_csi_pub']) * p['p_csi_ppm'],
        (1 - p['p_csi_pub']) * (1 - p['p_csi_ppm'])
    ])

    public = Sector(ent_pub, [alg1, alg2, alg3, alg4])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    txi = np.array([p['p_txi_pub'], p['p_txi_eng'], p['p_txi_pri']])

    return {
        'sys': sys,
        'p_txi': txi
    }


def get_intv_poc(p, target=0.8, p_txi_poc=0.95):
    sputum = Specimen('Sputum', p['p_loss_sputum'])
    ssm = Test('SSM', p['sens_ssm'], p['spec_ssm'], sputum)
    xpert = Test('Xpert', p['sens_xpert'], p['spec_xpert'], sputum)
    xpert_ssm = Test('Xpert_ss-', p['sens_xpert_ss-'], p['spec_xpert'], sputum)

    cdx = Test('CDx', p['sens_cdx'], p['spec_cdx'])

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_ssm, cdx=cdx)
    alg2 = Algorithm('SSM > CDx', ssm=ssm, cdx=cdx)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx)
    alg4 = Algorithm('CDx', cdx=cdx)

    ent_pub = np.array([
        p['p_ava_ssm_pub'] * p['p_ava_xpert_pub'],
        p['p_ava_ssm_pub'] * (1 - p['p_ava_xpert_pub']),
        (1 - p['p_ava_ssm_pub']) * p['p_ava_xpert_pub'],
        (1 - p['p_ava_ssm_pub']) * (1 - p['p_ava_xpert_pub'])
    ])

    p_sx, p_s, p_x, p_n = ent_pub
    p_sx, p_x = fillup(p_sx, p_x, target)
    p_s, p_x = fillup(p_s, p_x, target)
    p_n, p_x = fillup(p_n, p_x, target)
    ent_pub = np.array([p_sx, p_s, p_x, p_n])

    ent_eng = [p['p_ava_xpert_eng'], (1 - p['p_ava_xpert_eng'])]
    ent_eng[1], ent_eng[0] = fillup(ent_eng[1], ent_eng[0], target)

    ent = np.array([
        p['p_csi_pub'],
        (1 - p['p_csi_pub']) * p['p_csi_ppm'],
        (1 - p['p_csi_pub']) * (1 - p['p_csi_ppm'])
    ])

    public = Sector(ent_pub, [alg1, alg2, alg3, alg4])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    p_txi_pub = p['p_txi_pub']
    p_txi_eng = p['p_txi_eng']

    if p_txi_poc is not None:
        dx = np.array([alg.dx(1, 0).TruePos for alg in public.Algorithms])
        dx *= ent_pub
        dx /= dx.sum()
        p_txi_pub = (dx * np.array([p_txi_poc, p_txi_pub, p_txi_poc, p_txi_pub])).sum()

        dx = np.array([alg.dx(1, 0).TruePos for alg in engaged.Algorithms])
        dx *= ent_eng
        dx /= dx.sum()
        p_txi_eng = (dx * np.array([p_txi_poc, p_txi_eng])).sum()

    txi = np.array([p_txi_pub, p_txi_eng, p['p_txi_pri']])

    return {
        'sys': sys,
        'p_txi': txi
    }


def get_perfect_dx(p):
    cdx = Test('CDx', 1, 1)

    alg4 = Algorithm('CDx', cdx=cdx)

    ent = np.array([
        p['p_csi_pub'],
        (1 - p['p_csi_pub']) * p['p_csi_ppm'],
        (1 - p['p_csi_pub']) * (1 - p['p_csi_ppm'])
    ])

    public = Sector(np.array([1]), [alg4])
    engaged = Sector(np.array([1]), [alg4])
    private = Sector(np.array([1]), [alg4])

    sys = System(ent, public, engaged, private)

    p_txi_pub = p['p_txi_pub']
    p_txi_eng = p['p_txi_eng']

    txi = np.array([p_txi_pub, p_txi_eng, p['p_txi_pri']])

    return {
        'sys': sys,
        'p_txi': txi
    }


def get_perfect_txi(p):
    return {
        'sys': get_system(p),
        'p_txi': np.ones(3)
    }
