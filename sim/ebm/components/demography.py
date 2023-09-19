from sim.ebm.components.proc import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography(Process):
    def __init__(self, keys, inputs):
        Process.__init__(self, 'demo', keys=keys)
        self.Inputs = inputs
        self.Year0 = inputs.Year0

    def calculate_calc(self, t, y, pars, calc: dict, **kwargs):
        I = self.Keys

        dr_tb = np.zeros_like(y)

        r_die_tx_pub, r_die_tx_pri_u, r_die_tx_pri_i = pars['r_txd']

        # if 'intv' in kwargs and kwargs['intv'] is not None:
        #     r_die_tx_pub = kwargs['intv'].modify_td(t, r_die_tx_pub)

        dr_tb[I.Asym] = pars['r_die_asym']
        dr_tb[I.Sym] = pars['r_die_sym']
        dr_tb[I.ExCS] = pars['r_die_sym']
        dr_tb[I.TxPub] = r_die_tx_pub
        dr_tb[I.TxPriOnPub] = r_die_tx_pri_u
        dr_tb[I.TxPriOnPri] = r_die_tx_pri_i

        calc['die_tb'] = die_tb = dr_tb * y

        rates = self.Inputs(max(t, self.Year0))

        mu = rates['r_death'] - die_tb.sum(0) / y.sum(0)

        calc['die'] = mu * y
        calc['mig'] = rates['r_mig'] * y
        calc['bir'] = rates['r_birth'] * y.sum()
        if y.shape[1] > 1:
            calc['ageing'] = rates['r_ageing'] * y

    def compose_dya(self, ya, calc: dict):
        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        I = self.Keys

        dy -= calc['die_tb'] + calc['die'] + calc['mig']
        dy[I.U, 0] += calc['bir']

        da[I.A_Mor] += calc['die_tb'].sum()

        if y.shape[1] > 1:
            dy -= calc['ageing']
            dy[:, 1:] += calc['ageing'][:, :-1]

        return dy, da
