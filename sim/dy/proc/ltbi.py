import numpy as np
from sim.dy.proc.base import Process
from sim.dy.intervention.tx import get_intv_tx

__author__ = 'Chu-Chang Ku'
__all__ = ['LatentTB']


class LatentTB(Process):
    def __init__(self, keys):
        Process.__init__(self, 'ltbi', keys)
        self.BPaLM = get_intv_tx('BPaLM')

    def find_dya(self, t, y, da, pars, calc, intv=None):
        I = self.Keys
        dy = np.zeros_like(y)

        r_act, r_react = pars['r_act'], pars['r_react']

        r_rel, r_rel_te = pars['r_relapse'] * pars['k_relapse_adj'], pars['r_relapse_te']

        try:
            intv_vac = intv.Vac
            r_act, r_react = intv_vac.modify_prog(t, r_act, r_react)
            r_rel, r_rel_te = intv_vac.modify_rel(t, r_rel, r_rel_te)
        except AttributeError or KeyError:
            pass

        irr = pars['irr'] * np.exp(- pars['drt_act'] * max(t - pars['t0_decline'], 0))
        rr_rel_pub = pars['rr_relapse_pub'] / (0.3 + 0.7 * pars['rr_relapse_pub'])
        rr_rel_pri = 1 / (0.3 + 0.7 * pars['rr_relapse_pub'])

        r_rel_teu, r_rel_tei = r_rel_te * rr_rel_pub, r_rel_te * rr_rel_pri

        try:
            intv_tx = intv.Tx
            r_rel_teu, r_rel_tei = intv_tx.modify_rel(t, r_rel_teu, r_rel_tei)
        except AttributeError or KeyError:
            pass

        r_rel_teu, r_rel_tei = self.BPaLM.modify_rel(t, r_rel_teu, r_rel_tei)

        inc_act = irr * r_act * y[I.FLat]
        inc_react = irr * r_react * y[I.SLat]
        inc_rel_rhu = irr * r_rel_teu * y[I.RHighPub]
        inc_rel_stu = irr * r_rel * y[I.RStPub]
        inc_rel_rhi = irr * r_rel_tei * y[I.RHighPri]
        inc_rel_sti = irr * r_rel * y[I.RStPri]

        calc['inc_recent'] = inc_recent = inc_act
        calc['inc_retreated_u'] = inc_retreated_u = inc_rel_rhu + inc_rel_stu
        calc['inc_retreated_i'] = inc_retreated_i = inc_rel_rhi + inc_rel_sti
        calc['inc_remote'] = inc_remote = inc_react + inc_retreated_u + inc_retreated_i
        calc['inc'] = calc['inc_recent'] + calc['inc_remote']

        dy[I.FLat] -= inc_act
        dy[I.SLat] -= inc_react

        dy[I.RHighPub] -= inc_rel_rhu
        dy[I.RStPub] -= inc_rel_stu
        dy[I.RHighPri] -= inc_rel_rhi
        dy[I.RStPri] -= inc_rel_sti
        dy[I.Asym] = inc_recent + inc_remote

        da[I.A_IncRecent] = inc_recent.sum()
        da[I.A_IncRemote] = inc_remote.sum()
        da[I.A_Inc] = da[I.A_IncRecent] + da[I.A_IncRemote]
        da[I.A_IncTreatedPub] = inc_retreated_u.sum()
        da[I.A_IncTreated] = da[I.A_IncTreatedPub] + inc_retreated_i.sum()

        # Self-clearance / stabilisation
        r_clear = pars['r_clear']
        clear_sl = r_clear * y[I.SLat]
        clear_rstu = r_clear * y[I.RStPub]
        clear_rsti = r_clear * y[I.RStPri]

        r_lat, r_stab = pars['r_lat'], 1 / 1.5
        stab_fl = r_lat * y[I.FLat]
        stab_rhu = r_stab * y[I.RHighPub]
        stab_rhi = r_stab * y[I.RHighPri]

        dy[I.U] += clear_sl + clear_rstu + clear_rsti
        dy[I.FLat] -= stab_fl
        dy[I.SLat] += stab_fl - clear_sl
        dy[I.RHighPub] -= stab_rhu
        dy[I.RStPub] += stab_rhu - clear_rstu
        dy[I.RHighPri] -= stab_rhi
        dy[I.RStPri] += stab_rhi - clear_rsti
        return dy
