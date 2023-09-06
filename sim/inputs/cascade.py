import numpy as np
import numpy.random as rd
import json
from sim.healthcare.system import get_system, get_system_ts

__author__ = 'Chu-Chang Ku'
__all__ = ['CasRepo']


EXO = {
    'sens_ssm': 0.64,
    'spec_ssm': 0.98,
    'sens_xpert': 0.85,
    'sens_xpert_ss-': 0.64,
    'spec_xpert': 0.98,
    'dur_pub': 0.5,
    'p_loss_sputum': 0.15,
    'p_loss_swab': 0.02
}


class CasRepo:
    def __init__(self, pars, prev, txo):
        self.Pars = pars
        self.Prev = prev
        self.TxO = {
            'Succ': np.zeros(2),
            'Die': np.zeros(2),
        }

        for d in txo:
            if d['Index'] == 'TxSucc':
                if d['Tag'] == 'Pub':
                    self.TxO['Succ'][0] = d['M']
                else:
                    self.TxO['Succ'][1] = d['M']
            else:
                if d['Tag'] == 'Pub':
                    self.TxO['Die'][0] = d['M']
                else:
                    self.TxO['Die'][1] = d['M']

        self.TxO['LTFU'] = 1 - self.TxO['Succ'] - self.TxO['Die']

    def reform_pars(self, exo):
        i = 0
        while i <= 100:
            try:
                return self.prepare_pars(exo)
            except AssertionError:
                i += 1
        raise AssertionError('No valid parameter set')

    def prepare_pars(self, exo, pp=None):
        if pp is None:
            pp = rd.choice(self.Pars, 1)[0]

        p0 = dict(exo)
        p0.update(pp)
        pp = p0
        # pp = dict(pp)

        prev = self.Prev
        p_a, p_s, p_c = prev['PrevAsym'], prev['PrevSym'], prev['PrevExCS']

        mu_a = pp['r_die_asym'] + pp['r_sc']
        mu_s = mu_c = pp['r_die_sym'] + pp['r_sc']

        txis = np.array([pp[f'r_txi_{sec}'] for sec in ['pub', 'eng', 'pri']])
        txi = txis.sum()

        p_txi = np.array([pp[f'p_txi_{sec}'] for sec in ['pub', 'eng', 'pri']])

        det = txis / p_txi
        ppm = det[1] / det[1:].sum()
        p_ent_pub = pp['p_csi_pub']
        p_ent = np.array([p_ent_pub, (1 - p_ent_pub) * ppm, (1 - p_ent_pub) * (1 - ppm)])

        sys = pp['sys']
        cs = sys.seek_care(1, 0)
        p_dx = np.array([cs['Public'].TruePos, cs['Engaged'].TruePos, cs['Private'].TruePos])

        p_dx_all = (p_txi * p_dx * p_ent).sum()

        r_det = txi / p_c
        r_csi = p_c * (r_det + mu_c) / p_s
        r_onset = p_s * (r_csi + mu_s) / p_a

        txi0 = p_s * r_csi * p_dx_all

        assert txi > txi0
        txi1 = txi - txi0

        r_recsi = txi1 / p_c / p_dx_all

        p_dx0 = p_dx1 = p_dx

        alo = np.array([[1, 0], [pp['p_pri_on_pub'], 1 - pp['p_pri_on_pub']], [0, 1]])

        inc = (mu_a + r_onset) * p_a
        dur = np.array([0.5, pp['dur_pri']])

        ps = {
            'p_txi': p_txi,
            'p_dx': p_dx,
            'p_dx0': p_dx0,
            'p_dx1': p_dx1,
            'p_ent': p_ent,
            'r_onset': r_onset,
            'r_csi': r_csi,
            'r_recsi': r_recsi,
            'mu_a': mu_a,
            'mu_s': mu_s,
            'mu_c': mu_c,
            'tx_alo': alo,
            'tx_dur': dur,
            'r_txs': self.TxO['Succ'] / dur,
            'r_txd': self.TxO['Die'] / dur,
            'r_txl': self.TxO['LTFU'] / dur,
            'ppv': np.array([pp[f'ppv_{sec}'] for sec in ['pub', 'eng', 'pri']]),
            'inc': inc,
            'prev_ut': prev['PrevUt'],
            'prev_asc': (p_a, p_s, p_c)

        }
        ps.update(exo)

        if 'sys_ts' in pp:
            ps['sys_ts'] = pp['sys_ts']
            cs = sys.seek_care(1, 0)

        return ps

    # def reform(self, ps, rat_01=1):
    #     prev = self.Prev
    #     p_s, p_c = prev['PrevSym'], prev['PrevExCS']
    #
    #     ps1 = dict(ps)
    #
    #     r_csi, r_recsi, p_ent, p_txi = ps['r_csi'], ps['r_recsi'], ps['p_ent'], ps['p_txi']
    #
    #     pdx = (p_ent * p_txi * ps['p_dx']).sum()
    #
    #     txi = p_s * r_csi * pdx + p_c * r_recsi * pdx
    #
    #     k = (p_s * r_csi * rat_01 + p_c * r_recsi * 1) / txi
    #
    #     ps1['p_dx0'] = ps['p_dx'] / pdx * rat_01 / k
    #     ps1['p_dx1'] = ps['p_dx'] / pdx / k
    #     return ps1

    @staticmethod
    def load(file):
        with open(file, 'r') as f:
            js = json.load(f)

        pars = js['pars']
        ps = list()

        for p0 in pars:
            p = dict(EXO)
            p['sens_cdx'] = p0['sens_cdx']
            p['spec_cdx'] = p0['spec_cdx']
            p['p_ava_ssm_pub'] = p0['p_ava_ssm_pub']
            p['p_ava_xpert_pub'] = p0['p_ava_naat_pub']
            p['p_ava_xpert_eng'] = p0['p_ava_naat_eng']
            p['p_csi_ppm'] = p0['p_csi_ppm']
            p['p_csi_pub'] = p0['p_csi_pub']

            p['sys'] = get_system(p)

            if 'p_loss_swab' in p:
                p['sys_ts'] = get_system_ts(p)

            p.update(p0)
            ps.append(p)

        return CasRepo(ps, js['prev'], js['txo'])


if __name__ == '__main__':
    cr = CasRepo.load(
        '../../data/pars_cs_independent.json'
    )

    ps = cr.prepare_pars({
        'r_die_asym': 0.1,
        'r_die_sym': 0.2,
        'r_sc': 0.2
    })

    for k, v in ps.items():
        print(k, ': ', v)
