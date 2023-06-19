from sim.components.proc import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['ActiveTB']


class ActiveTB(Process):
    def __init__(self, keys):
        Process.__init__(self, 'atb', keys)

    def calculate_calc(self, t, y, pars, calc: dict, **kwargs):
        I = self.Keys

        calc['sc_a'] = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = pars['r_sc'] * y[I.ExCS]

        calc['onset'] = pars['r_onset'] * y[I.Asym]

        calc['txs'] = y[[I.TxPub, I.TxPri]] * pars['r_txs'].reshape((-1, 1))
        calc['txl'] = y[[I.TxPub, I.TxPri]] * pars['r_txl'].reshape((-1, 1))

        return calc

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        onset = calc['onset']
        sc_a, sc_s, sc_c = calc['sc_a'], calc['sc_s'], calc['sc_c']

        txs, txl = calc['txs'], calc['txl']

        dy[I.Asym] += - onset - sc_a
        dy[I.Sym] += onset - sc_s
        dy[I.ExCS] += - sc_c
        dy[I.SLat] += sc_a + sc_s + sc_c

        dy[I.TxPub] += - txs[0] - txl[0]
        dy[I.TxPri] += - txs[1] - txl[1]
        dy[I.RLow] += txs.sum(0)
        dy[I.RHigh] += txl.sum(0)

        return dy, da

    def measure(self, mea, t, y, pars, calc, **kwargs):
        pass
