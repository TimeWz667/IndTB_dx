from sim.dy.components.proc import Process
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
        calc['sc_r'] = pars['r_sc'] * y[I.ReCS]

        calc['onset'] = pars['r_onset'] * y[I.Asym]

        r_txo = 1 / pars['tx_dur']

        try:
            intv_acf = kwargs['intv'].ACF
            r_acf = intv_acf.modify_acf(t, 0)
        except AttributeError or KeyError:
            r_acf = 0

        calc['acf_a'] = r_acf * y[I.Asym]
        calc['acf_s'] = r_acf * y[I.Sym]
        calc['acf_c'] = r_acf * y[I.ExCS]
        calc['acf_r'] = r_acf * y[I.ReCS]

        calc['txo'] = y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri]] * r_txo.reshape((-1, 1, 1))
        calc['txso'] = y[[I.TxsPub, I.TxsPriOnPub]] * r_txo[:2].reshape((-1, 1, 1))

        r_acquire_dr = pars['r_acquire_dr'] if t > 2010 else 0

        calc['acquire_dr'] = y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], :, I.DS] * r_acquire_dr

        return calc

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        onset = calc['onset']
        sc_a, sc_s, sc_c, sc_r = calc['sc_a'], calc['sc_s'], calc['sc_c'], calc['sc_r']

        dy[I.Asym] += - onset - sc_a
        dy[I.Sym] += onset - sc_s
        dy[I.ExCS] += - sc_c
        dy[I.ReCS] += - sc_r
        dy[I.SLat] += sc_a + sc_s + sc_c + sc_r

        txo = calc['txo']
        dy[I.TxPub] -= txo[0]
        dy[I.TxPriOnPub] -= txo[1]
        dy[I.TxPriOnPri] -= txo[2]

        dy[I.RHighPub] += txo[0]
        dy[I.RHighPri] += txo[1] + txo[2]

        txo = calc['txso']
        dy[I.TxsPub] -= txo[0]
        dy[I.TxsPriOnPub] -= txo[1]

        dy[I.RHighPub] += txo[0]
        dy[I.RHighPri] += txo[1]


        acquire_dr = calc['acquire_dr']
        y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], :, I.DS] -= acquire_dr
        y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], :, I.DR] += acquire_dr

        acf_a, acf_s, acf_c, acf_r = calc['acf_a'], calc['acf_s'], calc['acf_c'], calc['acf_r']
        dy[I.Asym] -= acf_a
        dy[I.Sym] -= acf_s
        dy[I.ExCS] -= acf_c
        dy[I.ReCS] -= acf_r
        dy[I.TxPub] += acf_a + acf_s + acf_c + acf_r

        da[I.A_ACF] += (acf_a + acf_s + acf_c + acf_r).sum()

        return dy, da
