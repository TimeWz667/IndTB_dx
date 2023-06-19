import numpy as np
import sim.keys as I
from sim.ebm.dy import ModelPlain
from sim.inputs import load_inputs
from sims_pars.fit.targets import read_targets
__all__ = ['ModelRisk']


class ModelRisk(ModelPlain):
    def __init__(self, inputs):
        ModelPlain.__init__(self, 4, inputs)
        self.Prop = np.array([.7, .1, .1, .1])

    def reform_parameters(self, p):
        p = ModelPlain.reform_parameters(self, p)

        p['irr'] = irr = np.ones(4)

        irr[1] = 1 * p['irr_30']
        irr[2] = irr[1] * p['irr_20']
        irr[3] = irr[2] * p['irr_10']
        irr /= (irr * self.Prop).sum()

        return p

    def measure(self, t, ya, pars, intv=None):
        mea = ModelPlain.measure(self, t, ya, pars, intv=intv)
        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape((I.N_States, self.NDim))
        ns = y.sum(0)
        prev = (y[I.Asym] + y[I.Sym] + y[I.ExCS])

        for i in range(4):
            mea[f'N_T{i}'] = ns[i]
            mea[f'PrevUt_T{i}'] = prev[i] / ns[i]

        return mea

    def get_y0(self, pars):
        y0 = np.zeros((I.N_States, self.NDim))
        y0[I.Asym], y0[I.Sym], y0[I.ExCS] = pars['prev_asc']
        y0[I.SLat] = 0.5
        y0[I.U] = 1 - y0.sum(0)
        y0 *= self.Inputs.Demography.N0 * self.Prop.reshape((1, -1))
        a0 = np.zeros(I.N_Aux)
        return np.concatenate([y0.reshape(-1), a0])


if __name__ == '__main__':
    import matplotlib.pylab as plt
    import numpy.random as rd
    from sims_pars import bayes_net_from_script, sample

    rd.seed(1167)

    with open('../../data/prior.txt', 'r') as f:
        prior = bayes_net_from_script(f.read())

    exo0 = {
        'beta': 25,
        'rr_inf_asym': 0.8,
        'drt_trans': 0.05,
        'irr_30': 1.2,
        'irr_20': 1.2,
        'irr_10': 1.5
    }

    inp = load_inputs('../../data')
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    model0 = ModelRisk(inp)
    cr = model0.Inputs.Cascade
    ps = sample(prior, exo0)
    ps = cr.prepare_pars(ps)

    _, ms0, _ = model0.simulate_to_fit(ps, np.linspace(2000, 2030, 31))

    fig, axes = plt.subplots(2, 3)

    ms0.N.plot(ax=axes[0, 0])
    ms0.N_T0.plot(ax=axes[0, 0])
    ms0.N_T1.plot(ax=axes[0, 0])
    ms0.N_T2.plot(ax=axes[0, 0])
    ms0.N_T3.plot(ax=axes[0, 0])
    axes[0, 0].set_title('Population')

    ms0.PrevUt.plot(ax=axes[0, 1])
    ms0.PrevA.plot(ax=axes[0, 1])
    ms0.PrevS.plot(ax=axes[0, 1])
    ms0.PrevC.plot(ax=axes[0, 1])
    axes[0, 1].set_title('Prevalence')

    ms0.PrevUt_T0.plot(ax=axes[0, 2])
    ms0.PrevUt_T1.plot(ax=axes[0, 2])
    ms0.PrevUt_T2.plot(ax=axes[0, 2])
    ms0.PrevUt_T3.plot(ax=axes[0, 2])
    axes[0, 2].set_title('Prevalence, gp')

    ms0.NotiR.plot(ax=axes[1, 0])
    ms0.NotiPubR.plot(ax=axes[1, 0])
    ms0.NotiPriR.plot(ax=axes[1, 0])
    axes[1, 2].set_title('CNR')

    ms0.IncR.plot(ax=axes[1, 1])
    axes[1, 0].set_title('Incidence')
    ms0.MorR.plot(ax=axes[1, 2])
    axes[1, 1].set_title('Mortality')

    fig.tight_layout()
    plt.show()
