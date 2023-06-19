import numpy as np
from sim.components.proc import Process

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

        irr = pars['irr'] * np.exp(- pars['drt_act'] * max(t - pars['t0_decline'], 0))

        calc['inc_act'] = irr * r_act * y[I.FLat]
        calc['inc_react'] = irr * r_react * y[I.SLat]
        calc['inc_rel_rl'] = irr * r_rel_tc * y[I.RLow]
        calc['inc_rel_rh'] = irr * r_rel_td * y[I.RHigh]
        calc['inc_rel_st'] = irr * r_rel * y[I.RSt]

        calc['clear_sl'] = r_clear * y[I.SLat]
        calc['clear_rst'] = r_clear * y[I.RSt]

        calc['stab_fl'] = r_lat * y[I.FLat]
        calc['stab_rl'] = r_lat * y[I.RLow]
        calc['stab_rh'] = r_lat * y[I.RHigh]

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        dy[I.FLat] = - calc['inc_act'] - calc['stab_fl']
        dy[I.SLat] = calc['stab_fl'] - calc['inc_react'] - calc['clear_sl']
        dy[I.RLow] = - calc['inc_rel_rl'] - calc['stab_rl']
        dy[I.RHigh] = - calc['inc_rel_rh'] - calc['stab_rh']
        dy[I.RSt] = calc['stab_rl'] + calc['stab_rh'] - calc['inc_rel_st'] - calc['clear_rst']
        dy[I.Asym] = sum(calc[k] for k in ['inc_act', 'inc_react', 'inc_rel_rl', 'inc_rel_rh', 'inc_rel_st'])
        dy[I.U] = calc['clear_sl'] + calc['clear_rst']

        da[I.A_IncRecent] = sum(calc['inc_act'])
        da[I.A_IncRemote] = sum(calc['inc_react'] + calc['inc_rel_st'] + calc['inc_rel_rl'] + calc['inc_rel_rh'])
        da[I.A_Inc] = da[I.A_IncRecent] + da[I.A_IncRemote]

        return dy, da
