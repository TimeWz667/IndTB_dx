import numpy as np
import sim.keys as I
from sim.ebm.util import AbsModelODE
from sim.components import Demography, ActiveTB, LatentTB, Dx
from sim.inputs import load_inputs
__all__ = ['ModelPlain', 'ModelBaseline']


class ModelPlain(AbsModelODE):
    def __init__(self, n_dim, inputs):
        AbsModelODE.__init__(self, n_dim, inputs, 2022, 2030, dt=1, t_warmup=300, dfe=None)
        self.Year0 = inputs.Demography.Year0

        self.YearBaseline = 2022
        self.ProcDemo = Demography(I, inputs.Demography)
        self.ProcATB = ActiveTB(I)
        self.ProcDx = Dx(I)
        self.ProcLTBI = LatentTB(I)

    def reform_parameters(self, p):
        p = dict(p)
        p['sus'] = sus = np.zeros((I.N_States, self.NDim))
        sus[I.U] = 1
        sus[I.SLat] = p['rr_sus_slat']
        sus[I.RLow] = p['rr_sus_slat']
        sus[I.RHigh] = p['rr_sus_slat']
        sus[I.RSt] = p['rr_sus_slat']

        p['trans'] = trans = np.zeros((I.N_States, self.NDim))
        trans[I.Asym] = p['rr_inf_asym']
        trans[I.Sym] = 1
        trans[I.ExCS] = p['rr_inf_cs']
        p['irr'] = np.ones(self.NDim)

        return p

    def get_y0(self, pars):
        y0 = np.zeros((I.N_States, self.NDim))
        y0[I.Asym], y0[I.Sym], y0[I.ExCS] = pars['prev_asc']
        y0[I.SLat] = 0.5
        y0[I.U] = 1 - y0.sum(0)
        y0 *= self.Inputs.Demography.N0 / self.NDim
        a0 = np.zeros(I.N_Aux)
        return np.concatenate([y0.reshape(-1), a0])

    def calc_dy_transmission(self, t, y, pars, intv=None):
        foi = pars['beta'] * (y * pars['trans']).sum() / y.sum()

        if t > pars['t0_decline']:
            dt = max(self.Year0, max(t, pars['t0_decline'])) - self.Year0
            foi *= np.exp(- pars['drt_trans'] * dt)

        infection = pars['sus'] * foi * y
        dy = - infection
        dy[I.FLat] += infection.sum(0)
        return dy

    def __call__(self, t, ya, pars, intv=None):
        t = max(t, self.Year0)

        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape((I.N_States, self.NDim))

        dy, da = np.zeros_like(y), np.zeros_like(aux)

        dy += self.calc_dy_transmission(t, y, pars, intv=None)

        calc = dict()
        for proc in [self.ProcATB, self.ProcDx, self.ProcLTBI, self.ProcDemo]:
            (dy0, da0), calc = proc.find_dya(t, (y, aux), pars, calc=calc, intv=intv)
            dy += dy0
            da += da0

        if t <= self.Year0:
            dy -= dy.sum(0, keepdims=True) * y / y.sum(0, keepdims=True)

        return np.concatenate([dy.reshape(-1), da])

    def measure(self, t, ya, pars, intv=None):
        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape((I.N_States, self.NDim))
        n = y.sum()

        mea = dict(Year=t, N=n,
                   PrevUt=(y[I.Asym] + y[I.Sym] + y[I.ExCS]).sum() / n,
                   PrevA=y[I.Asym].sum() / n,
                   PrevS=y[I.Sym].sum() / n,
                   PrevC=y[I.ExCS].sum() / n,
                   PrevTxPub=y[I.TxPub].sum() / n,
                   PrevTxPri=y[I.TxPri].sum() / n,
                   LTBI=(y[I.LTBI].sum()) / n)

        mea['PrA'] = mea['PrevA'] / mea['PrevUt']
        mea['PrS'] = mea['PrevS'] / mea['PrevUt']
        mea['PrC'] = mea['PrevC'] / mea['PrevUt']

        mea['CumIncRecent'] = aux[I.A_IncRecent]
        mea['CumIncRemote'] = aux[I.A_IncRemote]
        mea['CumIncRetreat'] = aux[I.A_IncRetreat]
        mea['CumInc'] = aux[I.A_Inc]
        mea['CumMor'] = aux[I.A_Mor]
        mea['CumNotiPub'] = aux[I.A_NotiPub]
        mea['CumNotiPri'] = aux[I.A_NotiPri]
        mea['CumNoti'] = mea['CumNotiPub'] + mea['CumNotiPri']
        mea['CumYieldATB'] = aux[I.A_YieldATB]
        mea['CumYieldLTBI'] = aux[I.A_YieldLTBI]
        mea['CumYieldOT'] = aux[I.A_YieldOT]

        # calc = dict()
        # self.ProcDx.calculate_calc(t, y, pars, calc, intv=intv)
        #
        # mea['YieldATB'] = calc['acf_a'].sum() + calc['acf_s'].sum() + calc['acf_c'].sum()
        # mea['YieldLTBI'] = calc['acf_l'].sum() + calc['acf_td'].sum()
        # mea['YieldOT'] = calc['acf_ot'].sum()

        return mea


class ModelBaseline(ModelPlain):
    def __init__(self, inputs):
        ModelPlain.__init__(self, 1, inputs)


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
        'drt_trans': 0.02
    }

    inp = load_inputs('../../data')
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    model0 = ModelBaseline(inp)
    cr = model0.Inputs.Cascade
    ps = sample(prior, exo0)
    ps = cr.prepare_pars(ps)

    _, ms0, _ = model0.simulate_to_fit(ps, np.linspace(2000, 2030, 31))

    fig, axes = plt.subplots(2, 3)

    ms0.N.plot(ax=axes[0, 0])
    axes[0, 0].set_title('Population')

    ms0.PrevUt.plot(ax=axes[0, 1])
    ms0.PrevA.plot(ax=axes[0, 1])
    ms0.PrevS.plot(ax=axes[0, 1])
    ms0.PrevC.plot(ax=axes[0, 1])
    axes[0, 1].set_title('Prevalence')
    ms0.NotiR.plot(ax=axes[0, 2])
    ms0.NotiPubR.plot(ax=axes[0, 2])
    ms0.NotiPriR.plot(ax=axes[0, 2])
    axes[0, 2].set_title('CNR')

    ms0.IncR.plot(ax=axes[1, 0])
    axes[1, 0].set_title('Incidence')
    ms0.MorR.plot(ax=axes[1, 1])
    axes[1, 1].set_title('Mortality')

    fig.tight_layout()
    plt.show()
