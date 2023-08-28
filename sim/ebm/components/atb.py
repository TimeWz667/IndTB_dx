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

        if 'intv' in kwargs and kwargs['intv'] is not None:
            n_asym, n_sym = y[I.Asym], y[I.Sym] + y[I.ExCS]
            n = y.sum(0)
            r_acf_a, r_acf_s = kwargs['intv'].modify_acf(t, n_asym, n_sym, n)

            r_txs, r_txl = kwargs['intv'].modify_com(t, r_txs, r_txl)
        else:
            r_acf_a, r_acf_s = 0, 0

        calc['acf_a'] = r_acf_a * y[I.Asym]
        calc['acf_s'] = r_acf_s * y[I.Sym]
        calc['acf_c'] = r_acf_s * y[I.ExCS]

        calc['txs'] = y[[I.TxPub, I.TxPri]] * r_txs.reshape((-1, 1))
        calc['txl'] = y[[I.TxPub, I.TxPri]] * r_txl.reshape((-1, 1))

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
        dy[I.TxPri] += - txs[1] - txl[1]

        acf_a, acf_s, acf_c = calc['acf_a'], calc['acf_s'], calc['acf_c']
        dy[I.Asym] += - acf_a
        dy[I.Sym] += - acf_s
        dy[I.ExCS] += - acf_c

        dy[I.TxPub] += acf_a + acf_s + acf_c

        dy[I.RLow] += txs.sum(0)
        dy[I.RHigh] += txl.sum(0)

        dy[I.SLat] += sc_a + sc_s + sc_c

        da[I.A_ACF] += (acf_a + acf_s + acf_c).sum()

        return dy, da
