import numpy as np
import numpy.random as rd
import pickle as pkl
from sim.healthcare.system import get_system

__author__ = 'Chu-Chang Ku'
__all__ = ['CasRepo']


class CasRepo:
    def __init__(self, pars, prev):
        self.Pars = pars
        self.Prev = prev

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
            pp = rd.choice(self.Pars)

        pp = dict(pp)
        pp.update(exo)

        prev = self.Prev
        prev_a, prev_s, prev_c = prev['PrevAsym'], prev['PrevSym'], prev['PrevExCS']

        mu_a = pp['r_die_asym'] + pp['r_sc']
        mu_s = mu_c = pp['r_die_sym'] + pp['r_sc']

        p_ent, p_itt, p_dx, p_txi = pp['p_ent'], pp['p_itt'], pp['p_dx'], pp['p_txi']
        p_det = (p_ent * p_itt * p_dx * p_txi).sum()
        txi = pp['txi']

        r_det = txi / prev_c
        r_csi = (r_det + mu_c) * prev_c / prev_s
        r_onset = (r_csi + mu_s) * prev_s / prev_a
        inc = (r_onset + mu_a) * prev_a

        det0 = prev_s * r_csi * p_det
        det1 = txi - det0
        r_recsi = det1 / (prev_c * p_det)

        alo = np.array([
            [1, 0, 0],
            [0, pp['p_pri_on_pub'], 1 - pp['p_pri_on_pub']],
            [0, 0, 1]
        ])

        dur = np.array([0.5, 0.5, pp['dur_pri']])

        pp.update({
            'r_onset': r_onset,
            'r_csi': r_csi,
            'r_recsi': r_recsi,
            'r_det': r_det,
            'mu_a': mu_a,
            'mu_s': mu_s,
            'mu_c': mu_c,
            'tx_alo': alo,
            'tx_dur': dur,
            'inc': inc,
            'prev_ut': prev['PrevUt'],
            'prev_asc': (prev_a, prev_s, prev_c)
        })

        return pp

    @staticmethod
    def make_system(constructor, ps):
        return constructor(ps['src'])

    @staticmethod
    def load(file):
        with open(file, 'rb') as f:
            ps = pkl.load(f)

        return CasRepo(ps['Particles'], ps['Prev'])


if __name__ == '__main__':
    cr = CasRepo.load(
        '../../pars/pars_cas_cdx.pkl'
    )

    ps = cr.prepare_pars({
        'r_die_asym': 0.1,
        'r_die_sym': 0.2,
        'r_sc': 0.2
    })

    for k, v in ps.items():
        print(k, ': ', v)
