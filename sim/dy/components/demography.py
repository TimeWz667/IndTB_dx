from sim.dy.components.proc import Process
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

        dr_tb[I.Asym] = pars['r_die_asym']
        dr_tb[I.Sym] = pars['r_die_sym']
        dr_tb[I.ExCS] = pars['r_die_sym']
        dr_tb[I.ReCS] = pars['r_die_sym']

        for it in [I.TxPub, I.TxPriOnPub, I.TxPriOnPri, I.TxsPub, I.TxsPriOnPub]:
            dr_tb[it] = pars['r_die_tx']

        calc['die_tb'] = die_tb = dr_tb * y

        rates = self.Inputs(max(t, self.Year0))

        offset = die_tb.sum((0, 2)) / y.sum((0, 2))
        mu = rates['r_death'] - offset

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
