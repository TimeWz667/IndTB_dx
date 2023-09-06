from pydantic import BaseModel
from pydantic.types import confloat
from scipy.optimize import brentq
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class ACF(BaseModel):
    Coverage: confloat(ge=0, le=40) = 0
    Sens: confloat(ge=0, le=1) = 0.9


class DeRel(BaseModel):
    Target: confloat(ge=0, le=1) = 0


class Vac(BaseModel):
    Efficacy: confloat(ge=0, le=1) = 0
    Coverage: confloat(ge=0, le=0.9) = 0
    Year0: float = 2025
    Preflight: confloat(ge=0) = 2


class Swab(BaseModel):
    Uptake: confloat(ge=0, le=1) = 0
    XpertAccess: confloat(ge=0, le=1) = 0


class TxIni(BaseModel):
    Target: confloat(ge=0, le=1) = 0


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
    DeRel: DeRel = DeRel()
    Vac: Vac = Vac()
    Swab: Swab = Swab()
    TxIni: TxIni = TxIni()

    # Timeline for scaling up
    T0_Intv: float = 2025
    T1_Intv: float = 2027

    def get_wt(self, t):
        if t > self.T1_Intv:
            return 1
        elif t < self.T0_Intv:
            return 0
        else:
            return (t - self.T0_Intv) / (self.T1_Intv - self.T0_Intv)

    def modify_access(self, t, p_ent):
        return p_ent

    def modify_dx(self, t, p_dx0, p_dx1, sys_ts):
        if t > self.T0_Intv and self.Swab.Uptake > 0:
            if self.Swab.XpertAccess > 0:
                resid = self.Swab.XpertAccess * self.get_wt(t)
                ent = sys_ts.Public.Entry

                resid -= ent[0] + ent[2]
                if resid > 0:
                    ent1 = ent.copy()

                    if resid > ent[3]:
                        ent1[3] = 0
                        ent1[2] += ent[3]
                        ent1[1] -= (resid - ent[3])
                        ent1[0] += (resid - ent[3])
                    else:
                        ent1[3] -= resid
                        ent1[2] += resid
                    sys_ts.Public.EntryMask = ent1

                ent = sys_ts.Engaged.Entry

                self.Swab.XpertAccess * self.get_wt(t)

                if resid > ent[1]:
                    resid -= ent[1]
                    ent1 = ent.copy()

                    ent1[1] -= resid
                    ent1[0] += resid
                    sys_ts.Engaged.EntryMask = ent1

            cs = sys_ts.seek_care(1, 0)
            p_dx_ts = np.array([cs['Public'].TruePos, cs['Engaged'].TruePos, cs['Private'].TruePos])

            wt = self.get_wt(t) * self.Swab.Uptake

            p_dx0 = p_dx0 * (1 - wt) + p_dx_ts * wt
            p_dx1 = p_dx1 * (1 - wt) + p_dx_ts * wt

        return p_dx0, p_dx1

    def modify_acf(self, t, n_asym, n_sym, n):
        if t > self.T0_Intv and self.ACF.Coverage > 0:

            r_acf = self.ACF.Coverage * self.ACF.Sens
            return r_acf

        else:
            return 0

    def modify_cs(self, t, r_cs, r_rcs):
        # if t > self.T0_Intv and self.CS.Scale > 0:
        #     wt = self.get_wt(t)
        #
        #     r_cs1 = r_cs / max(1 - self.CS.Scale, 0.05)
        #     r_rcs1 = r_rcs / max(1 - self.CS.Scale, 0.05)
        #
        #     r_cs = r_cs + (r_cs1 - r_cs) * wt
        #     r_rcs = r_rcs + (r_rcs1 - r_rcs) * wt
        return r_cs, r_rcs

    def modify_td(self, t, r_evt):
        # if t > self.T0_Intv and self.TxDie.Scale > 0:
        #     wt = self.get_wt(t)
        #     r_evt1 = r_evt * (1 - self.TxDie.Scale)
        #     r_evt = r_evt + (r_evt1 - r_evt) * wt
        return r_evt

    def modify_txi(self, t, p_txi):
        if t > self.T0_Intv and self.TxIni.Target > 0:
            p = self.TxIni.Target
            p_txi1 = p_txi.copy()
            p_txi1[0] = max(p, p_txi[0])
            p_txi1[1] = max(p, p_txi[1])
            p_txi = p_txi + (p_txi1 - p_txi) * self.get_wt(t)
        return p_txi

    def modify_com(self, t, r_succ, r_ltfu):
        # if t > self.T0_Intv and self.TxCom.Scale > 0:
        #     wt = np.zeros_like(r_succ)
        #     wt[:2] = self.get_wt(t) * self.TxCom.Scale
        #     dif = r_ltfu * wt
        #     r_ltfu = r_ltfu - dif
        #     r_succ = r_succ + dif
        return r_succ, r_ltfu

    def modify_rel(self, t, r_rel, r_rel_tc, r_rel_tl, r_txs, r_txl):
        if t > self.T0_Intv and self.DeRel.Target > 0:
            r0 = (r_txs.sum() * r_rel_tc + r_txl.sum() * r_rel_tl) * 1.5 / (r_txs.sum() + r_txl.sum())

            wt = self.get_wt(t)

            rat = self.DeRel.Target / r0
            r_rel_tc1 = r_rel_tc * rat
            r_rel_tl1 = r_rel_tl * rat
            r_rel_tc = r_rel_tc + (r_rel_tc1 - r_rel_tc1) * wt
            r_rel_tl = r_rel_tl + (r_rel_tl1 - r_rel_tl) * wt
        return r_rel, r_rel_tc, r_rel_tl

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
        # if t > self.T0_Intv and self.TPT.Scale > 0:
        #     if t > self.T1_Intv + 2:
        #         apx_notified = (t - self.T1_Intv) * notif
        #     elif t < self.T1_Intv:
        #         apx_notified = (t - self.T0_Intv) * notif / 2
        #     else:
        #         apx_notified = ((t - self.T1_Intv) + (self.T1_Intv - self.T0_Intv) / 2) * notif
        #         apx_notified -= (t - 2 - self.T0_Intv) * notif / 2
        #
        #     contacts = apx_notified * 8 * 0.7 # 8 people household and 0.7 LTBI among contacts
        #     ipt = apx_notified * self.TPT.Scale * 0.7 * fl / (fl + sl)
        #     ipt = min(ipt, contacts)
        #
        #     k = r_act * fl / (fl + 7 * contacts)
        #     act0 = k * (fl - contacts) + 8 * k * (contacts - ipt) + 8 * k * (1 - 0.8) * ipt
        #     red = (act0 / fl) / r_act
        #     red = max(red, 0)
        #     r_act = r_act * red
        #     r_act_vac = r_act_vac * red

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
