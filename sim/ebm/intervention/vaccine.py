import numpy as np
import sim.ebm.keys as I

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_vac']


class IntvVac:
    def __init__(self, prot_inf, prot_prog, prot_rel, year0=2025, preflight=2):
        assert prot_rel >= 0
        assert prot_prog >= 0
        assert prot_inf >= 0
        self.Year0 = year0
        self.Preflight = preflight
        self.ProtInf = prot_inf
        self.ProtProg = prot_prog
        self.ProtRel = prot_rel

    def uptake(self, t):
        t0 = self.Year0
        t1 = t0 + self.Preflight

        if t > t1:
            return 1
        elif t < t0:
            return 0
        else:
            return (t - t0) / (t1 - t0)

    def modify_sus(self, t, sus):
        wt = self.uptake(t)
        if wt <= 0:
            return sus

        sus1 = sus.copy()

        if sus1.shape[1] > 1:
            k = (1 - wt) + wt * (1 - self.ProtInf)
            sus1[I.FLat, 2] *= k
            sus1[I.SLat, 2] *= k
            sus1[I.U, 2] *= k

        return sus1

    def modify_prog(self, t, r_act, r_react):
        wt = self.uptake(t)
        if wt <= 0:
            return r_act, r_react

        k = np.ones(8)
        k[2:5] = (1 - wt) + wt * (1 - self.ProtProg)

        return r_act * k, r_react * k

    def modify_rel(self, t, r_rel, r_rel_tc, r_rel_td):
        wt = self.uptake(t)
        if wt <= 0:
            return r_rel, r_rel_tc, r_rel_td

        k = (1 - wt) + wt * (1 - self.ProtRel)

        return r_rel * k, r_rel_tc * k, r_rel_td * k


def get_intv_vac(p, key):
    if key == 'BCG':
        return IntvVac(0.45, 0, 0, 2026)
    elif key == 'M72':
        return IntvVac(0, 0.6, 0, 2028)
    elif key == 'BCG-M72':
        return IntvVac(0.45, 0.6, 0, 2028)
    elif key == 'Recurrence':
        return IntvVac(0, 0, 0.5, 2028)

    return None

