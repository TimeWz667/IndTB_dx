import numpy as np
import sim.dy.keys as I
from sim.dy.util import AbsModelODE
from sim.dy.proc import Demography, Transmission, ActiveTB, LatentTB, Dx

__all__ = ['ModelBaseline']


class ModelBaseline(AbsModelODE):
    def __init__(self, inputs):
        n_agp = len(inputs.Demography.N0)
        AbsModelODE.__init__(self, (I.N_States, I.N_States_R, n_agp), inputs, 2023, 2041,
                             dt=1, t_warmup=300, dfe=None)
        self.Year0 = inputs.Demography.Year0
        self.N_Agp = n_agp
        self.YearBaseline = 2022
        self.ProcDemo = Demography(I, inputs.Demography)
        self.ProcTrans = Transmission(I)
        self.ProcATB = ActiveTB(I)
        self.ProcLTBI = LatentTB(I)
        self.ProcDx = Dx(I)

    def reform_parameters(self, p):
        p = dict(p)
        p['sus'] = sus = np.zeros(self.NDim)
        sus[I.U] = 1
        sus[[I.SLat, I.RHighPub, I.RStPub, I.RHighPri, I.RStPri]] = p['rr_sus_slat']

        p['trans_ds'] = trans = np.zeros(self.NDim)
        trans[I.Asym, I.DS] = p['rr_inf_asym']
        trans[[I.Sym, I.ExCS, I.ReCS], I.DS] = 1

        p['trans_dr'] = trans = np.zeros(self.NDim)
        trans[I.Asym, I.DR] = p['rr_inf_asym']
        trans[[I.Sym, I.ExCS, I.ReCS], I.DR] = 1

        p['trans_fr'] = trans = np.zeros(self.NDim)
        trans[I.Asym, I.FR] = p['rr_inf_asym']
        trans[[I.Sym, I.ExCS, I.ReCS], I.FR] = 1

        p['irr'] = np.ones(self.N_Agp)
        if self.N_Agp > 1:
            sus[:, :, 0] *= 0
            sus[:, :, 1] *= 0
            sus[:, :, 2] *= 1
            sus[:, :, 3] *= p['irr_25']
            sus[:, :, 4] *= p['irr_35']
            sus[:, :, 5] *= p['irr_45']
            sus[:, :, 6] *= p['irr_55']
            sus[:, :, 7] *= p['irr_65']

        return p

    def get_y0(self, pars):
        y0 = np.zeros(self.NDim)
        y0[I.Asym, 0], y0[I.Sym, :, 0], y0[I.ExCS, :, 0] = 0.0004, 0.001, 0.0001
        y0[I.SLat, 0] = 0.1
        y0[I.U, 0] = 1 - y0[:, 0].sum(0)
        y0 /= y0.sum((0, 1))
        y0 *= self.Inputs.Demography.N(self.Year0)
        a0 = np.zeros(I.N_Aux)
        return np.concatenate([y0.reshape(-1), a0])

    def __call__(self, t, ya, pars, intv=None):
        t = max(t, self.Year0)

        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape(self.NDim)

        dy, da = np.zeros_like(y), np.zeros_like(aux)

        calc = dict()
        for proc in [self.ProcDemo, self.ProcTrans, self.ProcATB, self.ProcLTBI, self.ProcDx]:
            dy0 = proc.find_dya(t, y, da, pars, calc=calc, intv=intv)
            dy += dy0

        if t <= self.Year0:
            dy -= dy.sum((0, 1), keepdims=True) * y / y.sum((0, 1), keepdims=True)
        return np.concatenate([dy.reshape(-1), da])

    def measure(self, t, ya, pars, intv=None):
        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape(self.NDim)
        n = y.sum()

        demo = self.Inputs.Demography(t)

        mea = dict(Year=t, N=n, N0=demo['N'].sum(),
                   PrevUt=(y[I.Asym] + y[I.Sym] + y[I.ExCS]).sum() / n, PrevA=y[I.Asym].sum() / n,
                   PrevS=y[I.Sym].sum() / n, PrevC=y[I.ExCS].sum() / n, PrevR=y[I.ReCS].sum() / n,
                   PrevTxPub=y[I.TxPub].sum() / n,
                   PrevTxPri=(y[I.TxPriOnPub].sum() + y[I.TxPriOnPri].sum()) / n, LTBI=(y[I.LTBI].sum()) / n)

        mea['PrA'] = mea['PrevA'] / mea['PrevUt']
        mea['PrS'] = mea['PrevS'] / mea['PrevUt']
        mea['PrC'] = mea['PrevC'] / mea['PrevUt']
        mea['PrevDR'] = y[I.UtTB, I.DR].sum() / y[I.UtTB].sum()

        mea['CumIncRecent'] = aux[I.A_IncRecent]
        mea['CumIncRemote'] = aux[I.A_IncRemote]
        mea['CumIncTreated'] = aux[I.A_IncTreated]
        mea['CumIncTreatedPub'] = aux[I.A_IncTreatedPub]

        mea['CumInc'] = aux[I.A_Inc]
        mea['CumMor'] = aux[I.A_Mor]
        mea['CumNotiPub'] = aux[I.A_NotiPub]
        mea['CumNotiPri'] = aux[I.A_NotiPri]
        mea['CumNoti'] = mea['CumNotiPub'] + mea['CumNotiPri']
        mea['CumACF'] = aux[I.A_ACF]

        return mea


if __name__ == '__main__':
    import matplotlib.pylab as plt
    import numpy.random as rd
    from sim.inputs import load_inputs
    from sims_pars import bayes_net_from_script, sample

    rd.seed(1166)

    exo0 = {
        'beta': 20,
        'rr_inf_asym': 0.8,
        'drt_trans': 0,
        'drt_act': 0,
        'rt_cs': 0.05,
        'rr_relapse_pub': 1,
        'k_relapse_adj': 1,
        'rr_beta_dr': 0.95,
        'rr_beta_fr': 0,
        'r_acquire_dr': 0.02,
        'irr_25': 1,
        'irr_35': 1,
        'irr_45': 1,
        'irr_55': 1,
        'irr_65': 1
    }

    with open('../../data/prior.txt', 'r') as f:
        scr = f.read()
    bn = bayes_net_from_script(scr)

    inp = load_inputs('../../pars', cs_suffix='cas_cdx', agp='who')
    inp.Demography.set_year0(2000)
    model0 = ModelBaseline(inp)

    cr = model0.Inputs.Cascade
    ps = sample(bn, cond=exo0)
    ps = cr.prepare_pars(ps)

    _, ms0, _ = model0.simulate_to_fit(ps, np.linspace(2005, 2030, 31))

    # print((ms0.IncTreatedPubR / ms0.IncTreatedR).tail(10))
    # print((ms0.IncTreatedR / ms0.IncR).tail(10))

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

    fig.tight_layout()
    plt.show()
