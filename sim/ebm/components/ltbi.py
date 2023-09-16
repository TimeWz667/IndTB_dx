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

            cov0 = (y[I.FLatVac] + y[I.SLatVac]).sum()
            cov0 = cov0 / (cov0 + (y[I.FLat] + y[I.SLat]).sum())

            r_vac, r_act_vac, r_react_vac = intv.modify_vac_act0(t, r_lat, r_act, r_react, cov0)
        else:
            r_vac, r_act_vac, r_react_vac = 0, r_act, r_react

        irr = pars['irr'] * np.exp(- pars['drt_act'] * max(t - pars['t0_decline'], 0))
        rr_rel_pub = pars['rr_relapse_pub'] / (0.3 + 0.7 * pars['rr_relapse_pub'])
        rr_rel_pri = 1 / (0.3 + 0.7 * pars['rr_relapse_pub'])

        calc['inc_act'] = irr * r_act * y[I.FLat]
        calc['inc_act_v'] = irr * r_act_vac * y[I.FLatVac]
        calc['inc_react'] = irr * r_react * y[I.SLat]
        calc['inc_react_v'] = irr * r_react * y[I.SLatVac]
        calc['inc_rel_rlu'] = irr * r_rel_tc * rr_rel_pub * y[I.RLowPub]
        calc['inc_rel_rhu'] = irr * r_rel_td * rr_rel_pub * y[I.RHighPub]
        calc['inc_rel_stu'] = irr * r_rel * rr_rel_pub * y[I.RStPub]
        calc['inc_rel_rli'] = irr * r_rel_tc * rr_rel_pri * y[I.RLowPri]
        calc['inc_rel_rhi'] = irr * r_rel_td * rr_rel_pri * y[I.RHighPri]
        calc['inc_rel_sti'] = irr * r_rel * rr_rel_pri * y[I.RStPri]

        calc['clear_sl'] = r_clear * y[I.SLat]
        calc['clear_slv'] = r_clear * y[I.SLatVac]
        calc['clear_rstu'] = r_clear * y[I.RStPub]
        calc['clear_rsti'] = r_clear * y[I.RStPri]

        calc['stab_fl'] = r_lat * y[I.FLat]
        calc['stab_flv'] = r_lat * y[I.FLatVac]
        calc['stab_rlu'] = r_lat * y[I.RLowPub]
        calc['stab_rhu'] = r_lat * y[I.RHighPub]
        calc['stab_rli'] = r_lat * y[I.RLowPri]
        calc['stab_rhi'] = r_lat * y[I.RHighPri]

        calc['vac_fl'] = r_vac * y[I.FLat]
        calc['vac_sl'] = r_vac * y[I.SLat]

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        inc_recent = calc['inc_act'] + calc['inc_act_v']
        inc_remote = calc['inc_react'] + calc['inc_react_v']
        inc_retreated_u = calc['inc_rel_stu'] + calc['inc_rel_rlu'] + calc['inc_rel_rhu']
        inc_retreated_i = calc['inc_rel_sti'] + calc['inc_rel_rli'] + calc['inc_rel_rhi']
        inc_remote += inc_retreated_u + inc_retreated_i

        dy[I.FLat] = - calc['inc_act'] - calc['stab_fl'] - calc['vac_fl']
        dy[I.SLat] = calc['stab_fl'] - calc['inc_react'] - calc['clear_sl'] - calc['vac_sl']

        dy[I.FLatVac] = calc['vac_fl'] - calc['inc_act_v'] - calc['stab_flv']
        dy[I.SLatVac] = calc['vac_sl'] + calc['stab_flv'] - calc['inc_react_v'] - calc['clear_slv']

        dy[I.RLowPub] = - calc['inc_rel_rlu'] - calc['stab_rlu']
        dy[I.RHighPub] = - calc['inc_rel_rhu'] - calc['stab_rhu']
        dy[I.RStPub] = calc['stab_rlu'] + calc['stab_rhu'] - calc['inc_rel_stu'] - calc['clear_rstu']
        dy[I.RLowPri] = - calc['inc_rel_rli'] - calc['stab_rli']
        dy[I.RHighPri] = - calc['inc_rel_rhi'] - calc['stab_rhi']
        dy[I.RStPri] = calc['stab_rli'] + calc['stab_rhi'] - calc['inc_rel_sti'] - calc['clear_rsti']
        dy[I.Asym] = inc_recent + inc_remote
        dy[I.U] = calc['clear_sl'] + calc['clear_slv'] + calc['clear_rstu'] + calc['clear_rsti']

        da[I.A_IncRecent] = inc_recent.sum()
        da[I.A_IncRemote] = inc_remote.sum()
        da[I.A_Inc] = da[I.A_IncRecent] + da[I.A_IncRemote]
        da[I.A_IncTreated] = (inc_retreated_u + inc_retreated_i).sum()
        da[I.A_IncTreatedPub] = inc_retreated_u.sum()
        return dy, da
