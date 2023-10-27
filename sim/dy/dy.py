import numpy as np
import sim.dy.keys as I
from sim.dy.util import AbsModelODE
from sim.dy.components import Demography, ActiveTB, LatentTB, Dx

__all__ = ['ModelPlain', 'ModelBaseline']


class ModelPlain(AbsModelODE):
    def __init__(self, n_dim, inputs):
        AbsModelODE.__init__(self, n_dim, inputs, 2023, 2041, dt=1, t_warmup=300, dfe=None)
        self.Year0 = inputs.Demography.Year0
        self.YearBaseline = 2022
        self.ProcDemo = Demography(I, inputs.Demography)
        self.ProcATB = ActiveTB(I)
        self.ProcDx = Dx(I)
        self.ProcLTBI = LatentTB(I)

    def reform_parameters(self, p):
        p = dict(p)
        p['sus'] = sus = np.zeros(self.NDim)
        sus[I.U] = 1
        sus[I.SLat] = p['rr_sus_slat']
        sus[I.RHighPub] = p['rr_sus_slat']
        sus[I.RStPub] = p['rr_sus_slat']
        sus[I.RHighPri] = p['rr_sus_slat']
        sus[I.RStPri] = p['rr_sus_slat']

        p['trans_ds'] = trans = np.zeros(self.NDim)
        trans[I.Asym, :, I.DS] = p['rr_inf_asym']
        trans[I.Sym, :, I.DS] = 1
        trans[I.ExCS, :, I.DS] = 1
        trans[I.ReCS, :, I.DS] = 1

        p['trans_dr'] = trans = np.zeros(self.NDim)
        trans[I.Asym, :, I.DR] = p['rr_inf_asym']
        trans[I.Sym, :, I.DR] = 1
        trans[I.ExCS, :, I.DR] = 1
        trans[I.ReCS, :, I.DR] = 1

        p['trans_fr'] = trans = np.zeros(self.NDim)
        trans[I.Asym, :, I.FR] = p['rr_inf_asym']
        trans[I.Sym, :, I.FR] = 1
        trans[I.ExCS, :, I.FR] = 1
        trans[I.ReCS, :, I.FR] = 1
        p['irr'] = np.ones(self.NDim[1:])

        return p

    def get_y0(self, pars):
        y0 = np.zeros(self.NDim)
        y0[I.Asym, :, 0], y0[I.Sym, :, 0], y0[I.ExCS, :, 0] = 0.0004, 0.001, 0.0001
        y0[I.SLat, :, 0] = 0.1
        y0[I.U, :, 0] = 1 - y0[:, :, 0].sum(0)
        y0 /= y0.sum()
        y0 *= self.Inputs.Demography.N0
        a0 = np.zeros(I.N_Aux)
        return np.concatenate([y0.reshape(-1), a0])

    def calc_dy_transmission(self, t, y, pars, intv=None):
        foi_s = pars['beta'] * (y * pars['trans_ds']).sum() / y.sum()
        foi_r = pars['beta'] * pars['rr_beta_dr'] * (y * pars['trans_dr']).sum() / y.sum()
        foi_f = pars['beta'] * pars['rr_beta_fr'] * (y * pars['trans_fr']).sum() / y.sum()

        if t > pars['t0_decline']:
            dt = max(self.Year0, max(t, pars['t0_decline'])) - self.Year0
            foi_s *= np.exp(- pars['drt_trans'] * dt)
            foi_r *= np.exp(- pars['drt_trans'] * dt)
            foi_f *= np.exp(- pars['drt_trans'] * dt)

        sus = pars['sus']
        if intv is not None:
            try:
                intv_vac = intv.Vac
                sus = intv_vac.modify_sus(t, sus)
            except AttributeError or KeyError:
                pass

        infection_s = sus * foi_s * y
        infection_r = sus * foi_r * y
        infection_f = sus * foi_f * y
        dy = - infection_s - infection_r - infection_f
        dy[I.FLat, :, I.DS] += infection_s.sum((0, 2))
        dy[I.FLat, :, I.DR] += infection_r.sum((0, 2))
        dy[I.FLat, :, I.FR] += infection_f.sum((0, 2))
        return dy

    def __call__(self, t, ya, pars, intv=None):
        t = max(t, self.Year0)

        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape(self.NDim)

        dy, da = np.zeros_like(y), np.zeros_like(aux)

        dy += self.calc_dy_transmission(t, y, pars, intv=intv)

        calc = dict()
        for proc in [self.ProcATB, self.ProcDx, self.ProcLTBI, self.ProcDemo]:
            (dy0, da0), calc = proc.find_dya(t, (y, aux), pars, calc=calc, intv=intv)
            dy += dy0
            da += da0

        if t <= self.Year0:
            dy -= dy.sum((0, 2), keepdims=True) * y / y.sum((0, 2), keepdims=True)
        return np.concatenate([dy.reshape(-1), da])

    def measure(self, t, ya, pars, intv=None):
        y, aux = ya[:-I.N_Aux], ya[-I.N_Aux:]
        y = y.reshape(self.NDim)
        n = y.sum()

        mea = dict(Year=t, N=n, PrevUt=(y[I.Asym] + y[I.Sym] + y[I.ExCS]).sum() / n, PrevA=y[I.Asym].sum() / n,
                   PrevS=y[I.Sym].sum() / n, PrevC=y[I.ExCS].sum() / n, PrevR=y[I.ReCS].sum() / n,
                   PrevTxPub=y[I.TxPub].sum() / n,
                   PrevTxPri=(y[I.TxPriOnPub].sum() + y[I.TxPriOnPri].sum()) / n, LTBI=(y[I.LTBI].sum()) / n)

        mea['PrA'] = mea['PrevA'] / mea['PrevUt']
        mea['PrS'] = mea['PrevS'] / mea['PrevUt']
        mea['PrC'] = mea['PrevC'] / mea['PrevUt']
        mea['PrevDR'] = y[I.UtTB, :, I.DR].sum() / y[I.UtTB].sum()

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


class ModelBaseline(ModelPlain):
    def __init__(self, inputs):
        ModelPlain.__init__(self, (I.N_States, 1, I.N_States_R), inputs)


if __name__ == '__main__':
    import matplotlib.pylab as plt
    import numpy.random as rd
    from sim.inputs import load_inputs
    from sims_pars import bayes_net_from_script, sample

    rd.seed(1166)

    exo0 = {
        'beta': 15,
        'rr_inf_asym': 0.8,
        'drt_trans': 0,
        'drt_act': 0,
        'rr_relapse_pub': 1,
        'k_relapse_adj': 1,
        'rr_beta_dr': 0.95,
        'rr_beta_fr': 0,
        'r_acquire_dr': 0.02
    }

    with open('../../data/prior.txt', 'r') as f:
        scr = f.read()
    bn = bayes_net_from_script(scr)

    inp = load_inputs('../../pars', cs_suffix='cas_cdx')
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    model0 = ModelBaseline(inp)

    cr = model0.Inputs.Cascade
    ps = sample(bn, cond=exo0)
    ps = cr.prepare_pars(ps)

    _, ms0, _ = model0.simulate_to_fit(ps, np.linspace(2000, 2030, 31))

    print((ms0.IncTreatedPubR / ms0.IncTreatedR).tail(10))
    print((ms0.IncTreatedR / ms0.IncR).tail(10))

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
