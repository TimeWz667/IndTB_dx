from sim.healthcare.diagnosis import *
import numpy as np
import numpy.random as rd
from functools import lru_cache

__author__ = 'Chu-Chang Ku'
__all__ = ['Sector', 'System', 'get_system']


class Sector:
    def __init__(self, p_ent, alg):
        self.Entry = p_ent
        self.EntryMask = None
        self.Algorithms = alg

    @lru_cache
    def seek_care(self, n_tb, n_nontb):
        res = Results()
        ent = self.Entry if self.EntryMask is None else self.EntryMask
        for i, p in enumerate(ent):
            res += self.Algorithms[i].dx(n_tb * p, n_nontb * p)
        return res

    def seek_care_sto(self, n_tb, n_nontb):
        ns_tb = rd.multinomial(n_tb, self.Entry)
        ns_nontb = rd.multinomial(n_nontb, self.Entry)
        res = Results()
        for i in range(len(self.Entry)):
            res += self.Algorithms[i].dx_sto(ns_tb[i], ns_nontb[i])
        return res


class System:
    def __init__(self, p_ent, pub, eng, pri):
        self.Entry = p_ent
        self.Public = pub
        self.Engaged = eng
        self.Private = pri

    def seek_care(self, n_tb, n_nontb):
        n_tb_pub, n_tb_eng, n_tb_pri = self.Entry * n_tb
        n_nontb_pub, n_nontb_eng, n_nontb_pri = self.Entry * n_nontb
        return {
            'Public': self.Public.seek_care(n_tb_pub, n_nontb_pub),
            'Engaged': self.Engaged.seek_care(n_tb_eng, n_nontb_eng),
            'Private': self.Private.seek_care(n_tb_pri, n_nontb_pri)
        }

    def seek_care_sto(self, n_tb, n_nontb):
        n_tb_pub, n_tb_eng, n_tb_pri = rd.multinomial(n_tb, self.Entry)
        n_nontb_pub, n_nontb_eng, n_nontb_pri = rd.multinomial(n_nontb, self.Entry)
        return {
            'Public': self.Public.seek_care_sto(n_tb_pub, n_nontb_pub),
            'Engaged': self.Engaged.seek_care_sto(n_tb_eng, n_nontb_eng),
            'Private': self.Private.seek_care_sto(n_tb_pri, n_nontb_pri)
        }


def get_system(p, pt):
    sputum = Specimen('Sputum', p['p_scanty'])
    ssm = Test('SSM', p['sens_s'], p['spec_s'], sputum)
    xpert = Test('Xpert', p['sens_x'], p['spec_x'], sputum)
    xpert_sn = Test('Xpert_ss-', p['sens_x_sn'], p['spec_x_sn'], sputum)
    cdx_bn = Test('CDx', p['sens_cdx_bn'], p['spec_cdx_bn'])

    cdx = Test('CDx', p['sens_cdx'], p['spec_cdx'])

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_sn, cdx=cdx_bn)
    alg2 = Algorithm('SSM > CDx', ssm=ssm, cdx=cdx_bn)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx_bn)
    alg4 = Algorithm('CDx', cdx=cdx)

    ent_pub = np.array([
        p['p_path_ssm'] * p['p_xpert_sn'],
        p['p_path_ssm'] * (1 - p['p_xpert_sn']),
        (1 - p['p_path_ssm'])
    ])

    ent_eng = np.array([
        p['p_xpert_eng'],
        (1 - p['p_xpert_eng'])
    ])

    txi = np.array([pt['txi_pub'], pt['txi_eng'], pt['txi_pri']])
    ptxi = np.array([0.95, 0.85, 0.85])
    pdx = np.array([p['pdx_pub'], p['pdx_eng'], p['pdx_pri']])

    p_ent = txi / (pdx * ptxi)
    p_ent /= p_ent.sum()

    public = Sector(ent_pub, [alg1, alg2, alg3, alg4])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    return System(p_ent, public, engaged, private)


if __name__ == '__main__':
    p0 = {
        'sens_s': 0.64,
        'spec_s': 0.98,
        'sens_x': 0.85,
        'sens_x_sn': 0.64,
        'spec_x': 0.98,
        'spec_x_sn': 0.64,
        'p_scanty': 0.15,
        'p_loss_swab': 0.05,
        'sens_cdx': 0.7,
        'spec_cdx': 0.95,
        'sens_cdx_bn': 0.7,
        'spec_cdx_bn': 0.95,
        'p_path_ssm': 0.8327,
        'p_xpert_sn': 0.9447,
        'p_xpert_eng': 0.1177,
        'pdx_pub': 0.9359,
        'pdx_eng': 0.7999,
        'pdx_pri': 0.7821,
    }

    system = get_system(p0, {'txi_pub': 1, 'txi_eng': 1, 'txi_pri': 1})

    system.Public.seek_care(1e4, 1e5).print()
    system.Public.seek_care_sto(1e4, 1e5).print()

    res = system.seek_care(1e4, 1e5)
    for k, v in res.items():
        print(k)
        v.print()

    print('---------------------------------------------')
    res = system.seek_care_sto(1e4, 1e5)
    for k, v in res.items():
        print(k)
        v.print()

    print('---------------------------------------------')
    print('Public')
    for alg in system.Public.Algorithms:
        print('Algorithm ', alg.Key)
        res0 = alg.dx(1000, 0)
        res1 = alg.dx(0, 1000)
        print('-- TP_Bac: ', res0.TruePos / 1000)
        print('-- FN_Bac: ', res0.FalseNeg / 1000)
        print('-- T_SSM: ', res0['N_test_SSM'] / 1000)
        print('-- T_Xpert: ', (res0['N_test_Xpert_ss-'] + res0['N_test_Xpert']) / 1000)
        print('-- FP_Bac: ', res1.FalsePos / 1000)
        print('-- F_SSM: ', res1['N_test_SSM'] / 1000)
        print('-- F_Xpert: ', (res1['N_test_Xpert_ss-'] + res1['N_test_Xpert']) / 1000)
