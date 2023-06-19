from scipy.interpolate import interp1d
from functools import lru_cache

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography:
    def __init__(self, src, year0=2000, mig=False):
        self.Source = src
        self.Years = src['Time']
        self.Year0 = year0

        self.YearRange = [min(self.Years), max(self.Years)]

        self.RateDeath = interp1d(self.Years, src['r_die'])
        self.RateBirth = interp1d(self.Years, src['r_bir'])
        self.RateMigration = interp1d(self.Years, src['r_mig'])
        self.N = interp1d(self.Years, src['Pop'])
        self.N0 = self.N(year0)
        self.HasMigration = mig

    def set_year0(self, year0):
        assert self.YearRange[0] < year0 < self.YearRange[1]
        self.Year0 = year0
        self.N0 = self.N(year0)

    @lru_cache(maxsize=1024)
    def __call__(self, time):
        if time < self.YearRange[0]:
            time = self.YearRange[0]
        elif time > self.YearRange[1]:
            time = self.YearRange[1]

        br, dr, mr = float(self.RateBirth(time)), float(self.RateDeath(time)), float(self.RateMigration(time))

        if not self.HasMigration:
            br, mr = br - mr, 0

        return {
            'Year': time,
            'r_birth': br,
            'r_die': dr,
            'r_mig': mr
        }


if __name__ == '__main__':
    import json

    with open('../../data/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src, mig=True)

    print('Year: 1990')
    print(demo(1990))

    print('Year: 2030')
    print(demo(2030))
