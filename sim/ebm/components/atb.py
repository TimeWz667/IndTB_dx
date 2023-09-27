from sim.ebm.components.proc import Process
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

        r_txs, r_txl = pars['r_txs'], pars['r_txl']

        # try:
        #     intv_tx = kwargs['intv'].Tx
        #     r_txs, r_txl = intv_tx.modify_txo(t, r_txs, r_txl)
        # except AttributeError or KeyError:
        #     pass

        try:
            intv_acf = kwargs['intv'].ACF
            r_acf = intv_acf.modify_acf(t, 0)
        except AttributeError or KeyError:
            r_acf = 0

        calc['acf_a'] = r_acf * y[I.Asym]
        calc['acf_s'] = r_acf * y[I.Sym]
        calc['acf_c'] = r_acf * y[I.ExCS]

        calc['txs'] = y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri]] * r_txs.reshape((-1, 1))
        calc['txl'] = y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri]] * r_txl.reshape((-1, 1))

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

        dy[I.TxPub] += - txs[0] - txl[0]
        dy[I.TxPriOnPub] += - txs[1] - txl[1]
        dy[I.TxPriOnPri] += - txs[2] - txl[2]

        acf_a, acf_s, acf_c = calc['acf_a'], calc['acf_s'], calc['acf_c']
        dy[I.Asym] += - acf_a
        dy[I.Sym] += - acf_s
        dy[I.ExCS] += - acf_c

        dy[I.TxPub] += acf_a + acf_s + acf_c

        # dy[I.RLowPub] += txs[0]
        # dy[I.RLowPri] += txs[1] + txs[2]
        dy[I.RHighPub] += txl[0] + txs[0]
        dy[I.RHighPri] += txl[1] + txl[2] + txs[1] + txs[2]

        dy[I.SLat] += sc_a + sc_s + sc_c

        da[I.A_ACF] += (acf_a + acf_s + acf_c).sum()

        return dy, da
