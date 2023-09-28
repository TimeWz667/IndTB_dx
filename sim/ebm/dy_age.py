from sim.ebm.dy import ModelPlain
import sim.ebm.keys as I

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelAgeGrp']


class ModelAgeGrp(ModelPlain):
    def __init__(self, inputs):
        ModelPlain.__init__(self, len(inputs.Demography.N0), inputs)

    def reform_parameters(self, p):
        p = ModelPlain.reform_parameters(self, p)
        sus = p['sus']
        sus[:, 0] *= 0
        sus[:, 1] *= 0
        sus[:, 2] *= p['irr_25']
        sus[:, 3] *= p['irr_35']
        sus[:, 4] *= p['irr_45']
        sus[:, 5] *= p['irr_55']
        sus[:, 6] *= p['irr_65']
        sus[:, 7] *= p['irr_65']
        return p

    def measure(self, t, ya, pars, intv=None):
        mea = ModelPlain.measure(self, t, ya, pars, intv=intv)

        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape((I.N_States, self.NDim))
        n = y.sum(0)

        calc = dict()
        self.ProcLTBI.calculate_calc(t, y, pars, calc, intv=intv)

        inc = calc['inc_act'] #+ calc['inc_act_v']
        inc += calc['inc_react'] #+ calc['inc_react_v']
        # inc += calc['inc_rel_stu'] + calc['inc_rel_rlu'] + calc['inc_rel_rhu']
        # inc += calc['inc_rel_sti'] + calc['inc_rel_rli'] + calc['inc_rel_rhi']
        inc += calc['inc_rel_stu'] + calc['inc_rel_rhu']
        inc += calc['inc_rel_sti'] + calc['inc_rel_rhi']

        ltbi = y[I.LTBI].sum(0)

        for i, lab in enumerate(self.Inputs.Demography.DimNames['Age']):
            mea[f'N_{lab}'] = n[i]
            mea[f'IncR_{lab}'] = inc[i] / n[i]
            mea[f'LTBI_{lab}'] = ltbi[i] / n[i]

        return mea


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pylab as plt
    import numpy.random as rd
    from sim.inputs import load_inputs
    from sims_pars import bayes_net_from_script, sample

    rd.seed(1166)

    exo0 = {
        'beta': 25,
        'rr_inf_asym': 0.8,
        'drt_trans': 0,
        'drt_act': 0,
        'rr_relapse_pub': 1.9,
        'irr_25': 2
    }

    with open('../../data/prior.txt', 'r') as f:
        scr = f.read()
    bn = bayes_net_from_script(scr)

    inp = load_inputs('../../pars', cs_suffix='bac_cdx_sector_2022_re', agp='who')
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    model0 = ModelAgeGrp(inp)

    cr = model0.Inputs.Cascade
    ps = sample(bn, cond=exo0)
    ps = cr.prepare_pars(ps)

    _, ms0, _ = model0.simulate_to_fit(ps, np.linspace(2000, 2030, 31))

    print(ms0.IncTreatedPubR / ms0.IncTreatedR)

    fig, axes = plt.subplots(2, 3)

    ms0.N.plot(ax=axes[0, 0])
    axes[0, 0].set_title('Population')

    ms0.PrevUt.plot(ax=axes[0, 1])
    ms0.PrevA.plot(ax=axes[0, 1])
    ms0.PrevS.plot(ax=axes[0, 1])
    ms0.PrevC.plot(ax=axes[0, 1])
    axes[0, 1].set_title('Prevalence')
    # ms0.NotiR.plot(ax=axes[0, 2])
    ms0.NotiPubR.plot(ax=axes[0, 2])
    ms0.NotiPriR.plot(ax=axes[0, 2])
    axes[0, 2].set_title('CNR')

    ms0.IncR.plot(ax=axes[1, 0])
    axes[1, 0].set_title('Incidence')
    ms0.MorR.plot(ax=axes[1, 1])
    axes[1, 1].set_title('Mortality')

    for agp in model0.Inputs.Demography.DimNames['Age']:
        ms0[f'IncR_{agp}'].plot(ax=axes[1, 2])
    axes[1, 2].set_title('Incidence, age')

    fig.tight_layout()
    plt.show()

