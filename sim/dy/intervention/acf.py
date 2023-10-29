import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['IntvACF', 'get_intv_acf']


class IntvACF:
    def __init__(self, pp, cov, sens, year0=2025, preflight=2):
        self.Coverage = cov
        self.Sensitivity = sens
        self.Year0 = year0
        self.Preflight = preflight

    def uptake(self, t):
        t0 = self.Year0
        t1 = t0 + self.Preflight

        if t > t1:
            return 1
        elif t < t0:
            return 0
        else:
            return (t - t0) / (t1 - t0)

    def modify_acf(self, t, r_acf):
        wt = self.uptake(t)
        if wt <= 0:
            return r_acf

        r_acf = wt * self.Coverage * self.Sensitivity
        return r_acf


def get_intv_acf(p, key):
    cov, sens = key.split("_")
    cov, sens = float(cov), float(sens)
    if cov <= 0:
        return None
    return IntvACF(p, cov, sens)
