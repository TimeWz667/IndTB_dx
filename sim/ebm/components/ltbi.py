import numpy as np
from sim.ebm.components.proc import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['LatentTB']


class LatentTB(Process):
    def __init__(self, keys):
        Process.__init__(self, 'ltbi', keys)

    def calculate_calc(self, t, y, pars, calc: dict, **kwargs):
        I = self.Keys

        r_clear = pars['r_clear']
        r_act, r_react = pars['r_act'], pars['r_react']
        r_rel, r_rel_tc, r_rel_td = pars['r_relapse'], pars['r_relapse_tc'], pars['r_relapse_td']
        r_lat = pars['r_lat']

        if 'intv' in kwargs and kwargs['intv'] is not None:
            intv = kwargs['intv']
            r_rel = intv.modify_rel(t, r_rel)

            cov0 = (y[I.FLatVac] + y[I.SLatVac]).sum()
            cov0 = cov0 / (cov0 + (y[I.FLat] + y[I.SLat]).sum())

            r_vac, r_act_vac, r_react_vac = intv.modify_vac_act0(t, r_lat, r_act, r_react, cov0)

            fl = y[[I.FLat, I.FLatVac]].sum()
            sl = y[[I.SLat, I.SLatVac]].sum()
            notif = (calc['tp0'] + calc['tp1'] + calc['fp'])[:2].sum()
            r_act, r_act_vac = intv.modify_tpt(t, r_act, r_act_vac, fl, sl, notif)
        else:
            r_vac, r_act_vac, r_react_vac = 0, r_act, r_react

        irr = pars['irr'] * np.exp(- pars['drt_act'] * max(t - pars['t0_decline'], 0))

        calc['inc_act'] = irr * r_act * y[I.FLat]
        calc['inc_act_v'] = irr * r_act_vac * y[I.FLatVac]
        calc['inc_react'] = irr * r_react * y[I.SLat]
        calc['inc_react_v'] = irr * r_react * y[I.SLatVac]
        calc['inc_rel_rl'] = irr * r_rel_tc * y[I.RLow]
        calc['inc_rel_rh'] = irr * r_rel_td * y[I.RHigh]
        calc['inc_rel_st'] = irr * r_rel * y[I.RSt]

        calc['clear_sl'] = r_clear * y[I.SLat]
        calc['clear_slv'] = r_clear * y[I.SLatVac]
        calc['clear_rst'] = r_clear * y[I.RSt]

        calc['stab_fl'] = r_lat * y[I.FLat]
        calc['stab_flv'] = r_lat * y[I.FLatVac]
        calc['stab_rl'] = r_lat * y[I.RLow]
        calc['stab_rh'] = r_lat * y[I.RHigh]

        calc['vac_fl'] = r_vac * y[I.FLat]
        calc['vac_sl'] = r_vac * y[I.SLat]

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        inc_recent = calc['inc_act'] + calc['inc_act_v']
        inc_remote = calc['inc_react'] + calc['inc_react_v']
        inc_remote += calc['inc_rel_st'] + calc['inc_rel_rl'] + calc['inc_rel_rh']

        dy[I.FLat] = - calc['inc_act'] - calc['stab_fl'] - calc['vac_fl']
        dy[I.SLat] = calc['stab_fl'] - calc['inc_react'] - calc['clear_sl'] - calc['vac_sl']

        dy[I.FLatVac] = calc['vac_fl'] - calc['inc_act_v'] - calc['stab_flv']
        dy[I.SLatVac] = calc['vac_sl'] + calc['stab_flv'] - calc['inc_react_v'] - calc['clear_slv']

        dy[I.RLow] = - calc['inc_rel_rl'] - calc['stab_rl']
        dy[I.RHigh] = - calc['inc_rel_rh'] - calc['stab_rh']
        dy[I.RSt] = calc['stab_rl'] + calc['stab_rh'] - calc['inc_rel_st'] - calc['clear_rst']
        dy[I.Asym] = inc_recent + inc_remote
        dy[I.U] = calc['clear_sl'] + calc['clear_slv'] + calc['clear_rst']

        da[I.A_IncRecent] = inc_recent.sum()
        da[I.A_IncRemote] = inc_remote.sum()
        da[I.A_Inc] = da[I.A_IncRecent] + da[I.A_IncRemote]
        return dy, da
