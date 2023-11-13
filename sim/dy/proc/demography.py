from sim.dy.proc.base import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography(Process):
    def __init__(self, keys, inp):
        Process.__init__(self, 'demo', keys=keys)
        self.Inputs = inp
        self.Year0 = inp.Year0
        try:
            self.N_Agp = len(inp.N0)
        except TypeError:
            self.N_Agp = 1

    def find_dya(self, t, y, da, pars, calc, intv=None):
        I = self.Keys

        dr_tb = np.zeros_like(y)

        dr_tb[I.Asym] = pars['r_die_asym']
        dr_tb[I.Sym] = pars['r_die_sym']
        dr_tb[I.ExCS] = pars['r_die_sym']
        dr_tb[I.ReCS] = pars['r_die_sym']

        dr_tb[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri]] = pars['r_die_tx']
        dr_tb[[I.TxNewPub, I.TxNewPriOnPub, I.TxNewPriOnPri]] = pars['r_die_tx']

        calc['die_tb'] = die_tb = dr_tb * y

        ns = y.sum((0, 1))
        rates = self.Inputs(max(t, self.Year0), ns)

        offset = die_tb.sum((0, 1)) / ns
        mu = rates['r_death'] - offset

        calc['die'] = mu * y
        calc['mig'] = rates['r_mig'] * y
        calc['bir'] = rates['r_birth'] * y.sum()
        if y.shape[1] > 1:
            calc['ageing'] = rates['r_ageing'] * y

        dy = np.zeros_like(y)
        dy -= calc['die_tb'] + calc['die'] - calc['mig']
        dy[I.U, 0, 0] += calc['bir']

        da[I.A_Mor] += calc['die_tb'].sum()

        if self.N_Agp > 1:
            dy -= calc['ageing']
            dy[:, :, 1:] += calc['ageing'][:, :, :-1]
        return dy
