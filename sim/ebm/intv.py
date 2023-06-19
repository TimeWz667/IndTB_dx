from pydantic import BaseModel
from pydantic.types import confloat, conint
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = []


class ACF(BaseModel):
    Yield: confloat(ge=0, le=40) = 0
    Asym: bool = False
    Target: str = 'All'  # All, 30%, 20%, 10%


class Intervention(BaseModel):
    ACF: ACF = ACF()

    # PrDx: PrDx = PrDx()
    T0_Intv: float = 2025
    T1_Intv: float = 2027

    def modify_acf(self, t, n_asym, n_sym, n):
        if t > self.T0_Intv and self.ACF.Yield > 0:

            n_tar = self.ACF.Yield * n.sum() * 1e-5

            if self.ACF.Asym:
                eli_a = eli_s = np.ones(4)
            else:
                eli_a = np.zeros(4)
                eli_s = np.ones(4)

            if self.ACF.Target == '10%':
                eli_s[:3] = 0
            elif self.ACF.Target == '20%':
                eli_s[:2] = 0
            elif self.ACF.Target == '30%':
                eli_s[:1] = 0

            r_acf = n_tar / (eli_a * n_asym + eli_s * n_sym).sum()
            r_acf = min(20, r_acf)

            return r_acf * eli_a, r_acf * eli_s

        else:
            return np.zeros_like(n_asym), np.zeros_like(n_sym)


if __name__ == '__main__':
    from sim.inputs import load_inputs
    from sim.ebm.dy_risk import ModelRisk
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

    inp = load_inputs('../../data')
    model0 = ModelRisk(inp)

    with open('../../data/prior.txt', 'r') as f:
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
