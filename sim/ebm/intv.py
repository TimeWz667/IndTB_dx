from pydantic import BaseModel
from pydantic.types import confloat
from typing import Union
from scipy.optimize import brentq
from sim.healthcare import System
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class ACF(BaseModel):
    Coverage: confloat(ge=0, le=40) = 0
    Sens: confloat(ge=0, le=1) = 0.9
    Year0: float = 2025
    Preflight: confloat(ge=0) = 2


class Tx(BaseModel):
    Year0: float = 2025
    Preflight: confloat(ge=0) = 2


class Vac(BaseModel):
    Efficacy: confloat(ge=0, le=1) = 0
    Coverage: confloat(ge=0, le=0.9) = 0
    Year0: float = 2025
    Preflight: confloat(ge=0) = 2


class Dx(BaseModel):
    System: str = "None"
    PrTxi: np.ndarray = np.zeros(3)
    Year0: float = 2025
    Preflight: confloat(ge=0) = 2

    class Config:
        arbitrary_types_allowed = True


def find_rate(cov, r_lat, r_act, r_react, r_act_v, r_react_v, dr):
    def fn(r):
        f0 = 1
        s0 = r_lat / (r + dr + r_react)
        fv0 = r / (r_lat + r_act_v + dr)
        sv0 = (r * s0 + r_lat * fv0) / (r_react_v + dr)
        return (fv0 + sv0) / (f0 + s0 + fv0 + sv0) - cov

    return brentq(fn, 0, 20)


class Intervention(BaseModel):
    ACF: ACF = ACF()
    Dx: Dx = Dx()
    Vac: Vac = Vac()
    Tx: Tx = Tx()

    def get_wt(self, t, sc):
        t0 = sc.Year0
        t1 = t0 + sc.Preflight

        if t > t1:
            return 1
        elif t < t0:
            return 0
        else:
            return (t - t0) / (t1 - t0)

    def modify_acf(self, t):
        if self.ACF.Year0 > 0 and self.ACF.Coverage > 0:
            r_acf = self.ACF.Coverage * self.ACF.Sens
            return r_acf
        else:
            return 0

    def modify_dx(self, t, p_ent, p_dx, p_txi, pars):
        if t > self.Dx.Year0 and self.Dx.System != 'None':
            wt = self.get_wt(t, self.Dx)
            sys = pars[self.Dx.System]['sys']
            test = sys.seek_care(1, 0)
            p_ent1 = sys.Entry
            p_dx1 = np.array([r.TruePos for r in test.values()]) / p_ent1
            p_ent = p_ent * (1 - wt) + p_ent1 * wt
            p_dx = p_dx * (1 - wt) + p_dx1 * wt
            p_txi = p_txi * (1 - wt) + self.Dx.PrTxi * wt

        return p_ent, p_dx, p_txi

    def modify_tx(self, t, r_txs, r_txl, r_txd):
        return r_txs, r_txl, r_txd

    def modify_vac_act0(self, t, r_lat, r_act, r_react, cov0):
        if t > self.Vac.Year0 and self.Vac.Coverage > 0:
            t0, t1 = self.Vac.Year0, self.Vac.Year0 + self.Vac.Preflight
            if t > t0:
                t = min(t, t1)
                cov = (t - t0) / (t1 - t0) * self.Vac.Coverage
                # cov = self.Vac.Coverage
                if cov > cov0 and cov > 0:
                    r_vac = 20 * (cov - cov0) / cov
                    r_vac = min(r_vac, r_lat * cov / (1 - cov))
                else:
                    r_vac = 0

                r_act = r_act * (1 - self.Vac.Efficacy)
                r_react = r_react * (1 - self.Vac.Efficacy)
            else:
                r_vac = 0
        else:
            r_vac = 0
        return r_vac, r_act, r_react


if __name__ == '__main__':
    from wr.inputs import load_inputs
    from wr.ebm.dy_risk import ModelRisk
    import matplotlib.pylab as plt
    import numpy.random as rd
    from sims_pars import bayes_net_from_script, sample

    rd.seed(1167)

    exo0 = {
        'beta': 25,
        'rr_inf_asym': 0.8,
        'drt_trans': 0.05,
        'irr_30': 1.2,
        'irr_20': 1.2,
        'irr_10': 1.5
    }

    inp = load_inputs('../../db_src/India')
    model0 = ModelRisk(inp)

    with open('../../db_src/prior.txt', 'r') as f:
        prior = bayes_net_from_script(f.read())

    cr = inp.Cascade
    ps = sample(prior, cond=exo0)
    ps = cr.prepare_pars(ps)

    y0, ms0, _ = model0.simulate_to_baseline(ps)

    y0 = y0.y[:, -1]

    _, ms_baseline, _ = model0.simulate_onward(y0, ps)

    intv = Intervention.parse_obj({
        'ACF': {
            'Yield': 10,
            'Asym': True,
            'Target': '20%'
        }
    })

    _, ms_intv1, _ = model0.simulate_onward(y0, ps, intv=intv)

    intv = Intervention.parse_obj({
        'ACF': {
            'Yield': 10,
            'Asym': False,
            'Target': '20%'
        }
    })

    _, ms_intv2, _ = model0.simulate_onward(y0, ps, intv=intv)

    fig, axes = plt.subplots(2, 3)

    ms_baseline.NotiR.plot(ax=axes[1, 0])
    ms_baseline.NotiPubR.plot(ax=axes[1, 0])
    ms_baseline.NotiPriR.plot(ax=axes[1, 0])
    ms_intv1.NotiR.plot(ax=axes[1, 0])
    ms_intv1.NotiPubR.plot(ax=axes[1, 0])
    ms_intv1.NotiPriR.plot(ax=axes[1, 0])
    ms_intv2.NotiR.plot(ax=axes[1, 0])
    ms_intv2.NotiPubR.plot(ax=axes[1, 0])
    ms_intv2.NotiPriR.plot(ax=axes[1, 0])
    axes[1, 0].set_title('CNR')

    ms_baseline.IncR.plot(ax=axes[1, 1])
    ms_intv1.IncR.plot(ax=axes[1, 1])
    ms_intv2.IncR.plot(ax=axes[1, 1])
    axes[1, 1].set_title('Incidence')
    ms_baseline.MorR.plot(ax=axes[1, 2])
    ms_intv1.MorR.plot(ax=axes[1, 2])
    ms_intv2.MorR.plot(ax=axes[1, 2])
    axes[1, 2].set_title('Mortality')

    fig.tight_layout()
    plt.show()
