import numpy as np
import numpy.random as rd
import json
from sim.healthcare.system import get_system

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
        ppm = pp['p_csi_ppm']
        p_ent_pub = pp['p_csi_pub']
        p_ent = np.array([p_ent_pub, (1 - p_ent_pub) * ppm, (1 - p_ent_pub) * (1 - ppm)])

        p_itt0 = np.array([pp[f'p_itt0_{sector}'] for sector in ['pub', 'eng', 'pri']])
        p_itt1 = np.array([pp[f'p_itt1_{sector}'] for sector in ['pub', 'eng', 'pri']])

        sys = pp['sys']
        cs = sys.seek_care(1, 0)
        p_dx = np.array([cs['Public'].TruePos, cs['Engaged'].TruePos, cs['Private'].TruePos])

        r_onset = pp['r_onset']
        r_csi = pp['r_csi']
        r_recsi = pp['r_recsi']

        alo = np.array([
            [1, 0, 0],
            [0, pp['p_pri_on_pub'], 1 - pp['p_pri_on_pub']],
            [0, 0, 1]
        ])

        inc = (mu_a + r_onset) * p_a
        dur = np.array([0.5, 0.5, pp['dur_pri']])

        ps = {
            'p_txi': p_txi,
            'p_dx': p_dx,
            'p_ent': p_ent,
            'p_itt0': p_itt0,
            'p_itt1': p_itt1,
            'r_onset': r_onset,
            'r_csi': r_csi,
            'r_recsi': r_recsi,
            'mu_a': mu_a,
            'mu_s': mu_s,
            'mu_c': mu_c,
            'tx_alo': alo,
            'tx_dur': dur,
            'r_txs': self.TxO['Succ'].repeat([2, 1]) / dur,
            'r_txd': self.TxO['Die'].repeat([2, 1]) / dur,
            'r_txl': self.TxO['LTFU'].repeat([2, 1]) / dur,
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
            p['p_ava_xpert_pub'] = p0['p_ava_xpert_pub']
            p['p_ava_xpert_eng'] = p0['p_ava_xpert_eng']
            p['p_csi_ppm'] = p0['p_csi_ppm']
            p['p_csi_pub'] = p0['p_csi_pub']

            p['sys'] = get_system(p)

            p.update(p0)
            ps.append(p)

        return CasRepo(ps, js['prev'], js['txo'])


if __name__ == '__main__':
    cr = CasRepo.load(
        '../../pars/pars_bac_cdx_sector_2022_re.json'
    )

    ps = cr.prepare_pars({
        'r_die_asym': 0.1,
        'r_die_sym': 0.2,
        'r_sc': 0.2
    })

    for k, v in ps.items():
        print(k, ': ', v)
