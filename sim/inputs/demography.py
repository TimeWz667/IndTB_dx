from scipy.interpolate import interp1d
from functools import lru_cache
import json

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography:
    def __init__(self, src, year0=2000):
        self.Source = src
        self.Years = src['Year']
        self.Year0 = year0

        self.YearRange = [min(self.Years), max(self.Years)]

        self.RateDeath = interp1d(self.Years, src['RateDeath'])
        self.RateBirth = interp1d(self.Years, src['RateBirth'])
        self.RateMigration = interp1d(self.Years, src['RateMig'])
        self.N = interp1d(self.Years, src['N'])
        self.N0 = self.N(year0)

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

        return {
            'Year': time,
            'r_birth': br,
            'r_die': dr,
            'r_mig': mr
        }

    @staticmethod
    def load(file, year0=2000):
        with open(file, 'r') as f:
            js = json.load(f)
        return Demography(js, year0=year0)


if __name__ == '__main__':
    demo = Demography.load('../../db_src/Delhi/pars_demo.json', year0=1990)

    print('Year: 1990')
    print(demo(1990))

    print('Year: 2030')
    print(demo(2030))
