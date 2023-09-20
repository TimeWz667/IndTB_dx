import numpy as np
from sim.healthcare.intv_dx import *

__author__ = 'Chu-Chang Ku'
__all__ = ['IntvDx', 'get_intv_dx']


class IntvDx:
    def __init__(self, pp, constructor, year0=2025, preflight=2, **kwargs):
        ss = constructor(pp['src'], **kwargs)
        self.System = ss['sys']
        self.PrTxi = ss['p_txi']
        self.Stats = self.System.seek_care(1, 0)
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

    def modify_dx(self, t, p_ent, p_dx, p_txi):
        wt = self.uptake(t)
        test = self.Stats

        p_ent1 = self.System.Entry
        p_dx1 = np.array([r.TruePos for r in test.values()]) / p_ent1
        p_ent = p_ent * (1 - wt) + p_ent1 * wt
        p_dx = p_dx * (1 - wt) + p_dx1 * wt
        p_txi = p_txi * (1 - wt) + self.PrTxi * wt
        return p_ent, p_dx, p_txi


def get_intv_dx(p, key, p_txi_poc=0.95):
    if key == 'TSwab':
        return IntvDx(p, get_intv_tswab)
    elif key == 'POC':
        return IntvDx(p, get_intv_poc, target=0.8, p_txi_poc=p_txi_poc)
    elif key == 'POC_Hi':
        return IntvDx(p, get_intv_poc, target=0.95, p_txi_poc=None)

    return None
