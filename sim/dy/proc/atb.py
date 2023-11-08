from sim.dy.proc.base import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['ActiveTB']


class ActiveTB(Process):
    def __init__(self, keys):
        Process.__init__(self, 'atb', keys)

    def find_dya(self, t, y, da, pars, calc, intv=None):
        I = self.Keys
        dy = np.zeros_like(y)

        sc_a = pars['r_sc'] * y[I.Asym]
        sc_s = pars['r_sc'] * y[I.Sym]
        sc_c = pars['r_sc'] * y[I.ExCS]
        sc_r = pars['r_sc'] * y[I.ReCS]
        calc['onset'] = onset = pars['r_onset'] * y[I.Asym]

        dy[I.Asym] += - onset - sc_a
        dy[I.Sym] += onset - sc_s
        dy[I.ExCS] += - sc_c
        dy[I.ReCS] += - sc_r
        dy[I.SLat] += sc_a + sc_s + sc_c + sc_r

        # Resistance acquiring
        r_acquire_dr = pars['r_acquire_dr'] if t > 2010 else 0

        acquire_dr = r_acquire_dr * y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.DS]
        dy[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.DS] -= acquire_dr
        dy[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.DR] += acquire_dr

        acquire_fr = r_acquire_dr * y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.DR]
        dy[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.DR] -= acquire_fr
        dy[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri], I.FR] += acquire_fr

        return dy
