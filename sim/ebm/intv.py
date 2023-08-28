from pydantic import BaseModel
from pydantic.types import confloat, conint
import numpy as np
from scipy.optimize import brentq

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class PPM(BaseModel):
    Scale: confloat(ge=0, le=1) = 0


class ImpDx(BaseModel):
    Dx: confloat(ge=0, le=1) = 0


class ACF(BaseModel):
    Yield: confloat(ge=0, le=40) = 0
    Asym: bool = False
    Target: str = 'All'  # All, 30%, 20%, 10%


class CS(BaseModel):
    Scale: confloat(ge=0.1, le=1) = 0


class DeRel(BaseModel):
    Scale: confloat(ge=0, le=1) = 0


class TPT(BaseModel):
    Scale: conint(ge=0, le=50) = 0


class Vac(BaseModel):
    Efficacy: confloat(ge=0, le=1) = 0
    Coverage: confloat(ge=0, le=0.9) = 0
    Year0: float = 2023
    Preflight: confloat(ge=0) = 2


class TxDie(BaseModel):
    Scale: confloat(ge=0, le=1) = 0


class TxCom(BaseModel):
    Scale: confloat(ge=0, le=1) = 0


def find_rate(cov, r_lat, r_act, r_react, r_act_v, r_react_v, dr):
    def fn(r):
        f0 = 1
        s0 = r_lat / (r + dr + r_react)
        fv0 = r / (r_lat + r_act_v + dr)
        sv0 = (r * s0 + r_lat * fv0) / (r_react_v + dr)

        return (fv0 + sv0) / (f0 + s0 + fv0 + sv0) - cov

    return brentq(fn, 0, 20)


class Intervention(BaseModel):
    PPM: PPM = PPM()
    ImpDx: ImpDx = ImpDx()
    ACF: ACF = ACF()
    CS: CS = CS()
    TPT: TPT = TPT()
    DeRel: DeRel = DeRel()
    Vac: Vac = Vac()
    TxDie: TxDie = TxDie()
    TxCom: TxCom = TxCom()
    # PrDx: PrDx = PrDx()

    # Timeline for scaling up
    T0_Intv: float = 2023
    T1_Intv: float = 2025

    def get_wt(self, t):
        if t > self.T1_Intv:
            return 1
        elif t <= self.T0_Intv:
            return 0
        else:
            return (t - self.T0_Intv) / (self.T1_Intv - self.T0_Intv)

    def modify_access(self, t, p_ent):
        if t > self.T0_Intv and self.PPM.Scale > 0:
            wt = self.get_wt(t)
            p_ent1 = p_ent.copy()

            temp = p_ent[2] * self.PPM.Scale
            p_ent1[1] += temp
            p_ent1[2] -= temp

            p_ent = p_ent + wt * (p_ent1 - p_ent)
        return p_ent

    def modify_dx(self, t, p_dx0, p_dx1):
        if t > self.T0_Intv and self.ImpDx.Dx > 0:
            wt = np.zeros(3)
            wt[:2] = self.get_wt(t)

            p_dx0 = p_dx0 + (1 - p_dx0) * wt * self.ImpDx.Dx
            p_dx1 = p_dx1 + (1 - p_dx1) * wt * self.ImpDx.Dx

        return p_dx0, p_dx1

    def modify_acf(self, t, n_asym, n_sym, n):
        if t > self.T0_Intv and self.ACF.Yield > 0:

            n_tar = self.ACF.Yield * n.sum() * 1e-5

            if self.ACF.Asym:
                eli_a = np.ones(4)
                eli_s = np.ones(4)
            else:
                eli_a = np.zeros(4)
                eli_s = np.ones(4)

            if self.ACF.Target == '10%':
                eli_a[:3] = eli_s[:3] = 0
            elif self.ACF.Target == '20%':
                eli_a[:2] = eli_s[:2] = 0
            elif self.ACF.Target == '30%':
                eli_a[:1] = eli_s[:1] = 0

            r_acf = n_tar / (eli_a * n_asym + eli_s * n_sym).sum()
            r_acf = min(20, r_acf)

            return r_acf * eli_a, r_acf * eli_s

        else:
            return np.zeros_like(n_asym), np.zeros_like(n_sym)

    def modify_cs(self, t, r_cs, r_rcs):
        if t > self.T0_Intv and self.CS.Scale > 0:
            wt = self.get_wt(t)

            r_cs1 = r_cs / max(1 - self.CS.Scale, 0.05)
            r_rcs1 = r_rcs / max(1 - self.CS.Scale, 0.05)

            r_cs = r_cs + (r_cs1 - r_cs) * wt
            r_rcs = r_rcs + (r_rcs1 - r_rcs) * wt
        return r_cs, r_rcs

    def modify_td(self, t, r_evt):
        if t > self.T0_Intv and self.TxDie.Scale > 0:
            wt = self.get_wt(t)
            r_evt1 = r_evt * (1 - self.TxDie.Scale)
            r_evt = r_evt + (r_evt1 - r_evt) * wt
        return r_evt

    def modify_com(self, t, r_succ, r_ltfu):
        if t > self.T0_Intv and self.TxCom.Scale > 0:
            wt = np.zeros_like(r_succ)
            wt[:2] = self.get_wt(t) * self.TxCom.Scale
            dif = r_ltfu * wt
            r_ltfu = r_ltfu - dif
            r_succ = r_succ + dif
        return r_succ, r_ltfu

    def modify_rel(self, t, r_rel):
        if t > self.T0_Intv and self.DeRel.Scale > 0:
            wt = self.get_wt(t)
            r_rel1 = (1 - self.DeRel.Scale) * r_rel
            r_rel = r_rel + (r_rel1 - r_rel) * wt
        return r_rel

    def modify_vac_act0(self, t, r_lat, r_act, r_react, cov0):
        if t > self.T0_Intv and self.Vac.Coverage > 0:
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

    # def modify_vac_act(self, t, r_lat, r_act, r_react, r_die):
    #     r_act_v, r_react_v = r_act, r_react
    #     if t > self.T0_Intv and self.Vac.Coverage > 0:
    #         t0, t1 = self.Vac.Year0, self.Vac.Year0 + self.Vac.Preflight
    #         if t > t0:
    #             t = min(t, t1)
    #             cov = self.Vac.Coverage
    #
    #             r_act_v = r_act * (1 - self.Vac.Efficacy)
    #             r_react_v = r_react * (1 - self.Vac.Efficacy)
    #             try:
    #                 r_vac = find_rate(cov, r_lat, r_act[0], r_react[0], r_act_v[0], r_react_v[0], r_die)
    #             except ValueError:
    #                 r_vac = 0
    #         else:
    #             r_vac = 0
    #     else:
    #         r_vac = 0
    #     return r_vac * 20, r_act_v, r_react_v

    def modify_tpt(self, t, r_act, r_act_vac, fl, sl, notif):
        if t > self.T0_Intv and self.TPT.Scale > 0:
            if t > self.T1_Intv + 2:
                apx_notified = (t - self.T1_Intv) * notif
            elif t < self.T1_Intv:
                apx_notified = (t - self.T0_Intv) * notif / 2
            else:
                apx_notified = ((t - self.T1_Intv) + (self.T1_Intv - self.T0_Intv) / 2) * notif
                apx_notified -= (t - 2 - self.T0_Intv) * notif / 2

            contacts = apx_notified * 8 * 0.7 # 8 people household and 0.7 LTBI among contacts
            ipt = apx_notified * self.TPT.Scale * 0.7 * fl / (fl + sl)
            ipt = min(ipt, contacts)

            k = r_act * fl / (fl + 7 * contacts)
            act0 = k * (fl - contacts) + 8 * k * (contacts - ipt) + 8 * k * (1 - 0.8) * ipt
            red = (act0 / fl) / r_act
            red = max(red, 0)
            r_act = r_act * red
            r_act_vac = r_act_vac * red

        return r_act, r_act_vac


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
