from scipy.interpolate import interp1d
import numpy as np


__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


def interp_cont(years, v):
    return interp1d(years, v, axis=0, kind='linear', bounds_error=False, fill_value=(v[0], v[-1]))


class Demography:
    def __init__(self, src):
        self.DimNames = src['dimnames']
        self.Source = dict(src)
        self.Years = years = src['Year']
        self.YearSpan = min(years), max(years)
        self.Year0 = self.YearSpan[0] + 0.5

        fn = interp_cont

        self.N = fn(years, src['N'])
        self.RateBirth = fn(years, src['RateBirth'])
        self.RateDeath = fn(years, src['RateDeath'])
        self.RateAgeing = fn(years, src['RateAgeing'])
        self.RateMig = fn(years, src['RateMig']) if 'RateMig' in src else None

        self.N0 = self.N(self.Year0)

    def set_year0(self, year0):
        assert self.YearSpan[0] < year0 < self.YearSpan[1]
        self.Year0 = year0
        self.N0 = self.N(year0)

    def calc_mig(self, t, y):
        n = self.N(t)
        mr = 50 * (n - y) / n
        return mr

    def __call__(self, t, y=None):
        t = max(t, self.Year0)

        pars = {
            'N': self.N(t),
            'r_birth': self.RateBirth(t),
            'r_ageing': self.RateAgeing(t),
            'r_death': self.RateDeath(t)
        }
        if self.RateMig is not None:
            pars['r_mig'] = self.RateMig(t)
        elif y is not None:
            pars['r_mig'] = self.calc_mig(t, y)
        else:
            pars['r_mig'] = np.zeros_like(pars['r_death'])

        return pars


if __name__ == '__main__':
    import pickle as pkl

    with open('../../pars/ind_who_70to35.pkl', 'rb') as f:
        src = pkl.load(f)

    demo = Demography(src)
    for k, v in demo(2020).items():
        print(k, v)

    print('----------------------------------------------------')
    with open('../../pars/ind_all_70to35.pkl', 'rb') as f:
        src = pkl.load(f)

    demo = Demography(src)
    for k, v in demo(2020).items():
        print(k, v)
