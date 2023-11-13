import numpy as np
from sim.healthcare.intv_dx import *

__author__ = 'Chu-Chang Ku'
__all__ = ['IntvPPM', 'get_intv_ppm', 'IntvDx', 'get_intv_dx']


class IntvPPM:
    def __init__(self, target, year0=2025, preflight=2):
        self.Target = target
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

    def modify_ent(self, t, p_ent):
        wt = self.uptake(t)

        if wt <= 0:
            return p_ent

        p_ent1 = p_ent.copy()
        ppm0 = p_ent1[1] / (p_ent1[1] + p_ent1[2])

        if ppm0 > self.Target:
            return p_ent

        ppm1 = (1 - wt) * ppm0 + wt * self.Target
        p_ent1[1], p_ent1[2] = ppm1 * (p_ent1[1] + p_ent1[2]), (1 - ppm1) * (p_ent1[1] + p_ent1[2])
        return p_ent1


def get_intv_ppm(target):
    if target >= 0:
        return IntvPPM(target)
    return None


class IntvDx:
    def __init__(self, pp, constructor, year0=2025, preflight=2, redgap_itt=0.0, **kwargs):
        ss = constructor(pp, **kwargs)
        self.System = ss['sys']
        self.PrTxi = ss['p_txi']
        self.RedGapItt = redgap_itt

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

    def modify_dx(self, t, p_itt, p_ent, p_dx, p_txi):
        wt = self.uptake(t)
        test = self.Stats

        if self.RedGapItt > 0:
            p_itt1 = np.ones_like(p_ent) * p_itt
            p_itt1[:2] = 1 - (1 - p_itt1[:2]) * (1 - self.RedGapItt)
            p_itt = p_itt * (1 - wt) + p_itt1 * wt

        p_ent1 = self.System.Entry

        p_dx1 = np.array([r.TruePos for r in test.values()]) / p_ent1
        p_ent = p_ent * (1 - wt) + p_ent1 * wt
        p_dx = p_dx * (1 - wt) + p_dx1 * wt
        p_txi = p_txi * (1 - wt) + self.PrTxi * wt
        return p_itt, p_ent, p_dx, p_txi


def get_intv_dx(p, key, p_txi_poc=0.95):
    if key == 'TSwab':
        return IntvDx(p, get_intv_tswab)
    elif key == 'TSwab_ITT':
        return IntvDx(p, get_intv_tswab, redgap_itt=0.5)
    elif key == 'POC':
        return IntvDx(p, get_intv_poc, target=0.8, p_txi_poc=p_txi_poc)
    elif key == 'POC_ITT':
        return IntvDx(p, get_intv_poc, target=0.8, p_txi_poc=p_txi_poc, redgap_itt=0.5)
    elif key == 'POC_Hi':
        return IntvDx(p, get_intv_poc, target=0.95, p_txi_poc=None)
    elif key == 'POC_Hi_ITT':
        return IntvDx(p, get_intv_poc, target=0.95, p_txi_poc=None, redgap_itt=0.5)
    elif key == 'PerfectDx':
        return IntvDx(p, get_perfect_dx, redgap_itt=0.5)
    elif key == 'PerfectTxi':
        return IntvDx(p, get_perfect_txi, redgap_itt=0.5)

    return None
