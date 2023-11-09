import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['IntvACF', 'get_intv_acf']


class IntvACF:
    def __init__(self, cov, sens, new_tx=False, year0=2025, preflight=2):
        self.Coverage = cov
        self.Sensitivity = sens
        self.UseNewTx = new_tx
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

    def modify_acf(self, t, r_acf, r_acf_new):
        wt = self.uptake(t)
        if wt <= 0:
            return r_acf, r_acf_new

        r_acf = wt * self.Coverage * self.Sensitivity

        if self.UseNewTx:
            return 0, r_acf
        else:
            return r_acf, 0


def get_intv_acf(key):
    cov, sens, new = key.split("_")
    cov, sens = float(cov), float(sens)
    new = new == '1'
    if cov <= 0:
        return None
    return IntvACF(cov, sens, new)
