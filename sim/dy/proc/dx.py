from sim.dy.proc.base import Process
import numpy as np
from sim.dy.intervention.tx import get_intv_tx

__author__ = 'Chu-Chang Ku'
__all__ = ['Dx']


class Dx(Process):
    def __init__(self, keys):
        Process.__init__(self, 'dx', keys)
        self.BPalM = get_intv_tx('BPaLM')

    def find_dya(self, t, y, da, pars, calc, intv=None):
        I = self.Keys
        dy = np.zeros_like(y)

        k = np.exp(pars['rt_cs'] * max(t - pars['t0_decline'], 0))
        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']
        r_csi *= k
        r_recsi *= k

        p_ent, p_itt, p_pdx, p_txi = pars['p_ent'], pars['p_itt'], pars['p_dx'], pars['p_txi']
        p_dx_xpert = pars['p_dx_xpert']

        try:
            p_itt, p_ent, p_pdx, p_txi = intv.Dx.modify_dx(t, p_itt, p_ent, p_pdx, p_txi)
        except AttributeError or KeyError:
            pass

        try:
            p_ent = intv.PPM.modify_ent(t, p_ent)
        except AttributeError or KeyError:
            pass

        p_txi = self.BPalM.modify_txi(t, p_txi)
        try:
            p_txi = intv.Tx.modify_txi(t, p_txi)
        except AttributeError or KeyError:
            pass

        #
        pdx = p_ent * p_itt * p_pdx * p_txi
        p_fn = 1 - pdx.sum()
        ppv, alo = pars['ppv'], pars['tx_alo']
        p_cure0, p_cure0_dr = pars['p_cure'], pars['p_cure_dr']
        p_cure0_u, p_cure0_e, p_cure0_i = p_cure0, p_cure0, p_cure0
        p_cure0_dr_u, p_cure0_dr_e, p_cure0_dr_i = p_cure0_dr, p_cure0_dr, p_cure0_dr

        p_cure0_u, p_cure0_e, p_cure0_dr_u, p_cure0_dr_e = self.BPalM.modify_cure(t, p_cure0_u, p_cure0_e,
                                                                                     p_cure0_dr_u, p_cure0_dr_e)
        try:
            intv_tx = intv.Tx
            p_cure0_u, p_cure0_e, p_cure0_dr_u, p_cure0_dr_e = intv.Tx.modify_cure(t, p_cure0_u, p_cure0_e,
                                                                                   p_cure0_dr_u, p_cure0_dr_e)
        except AttributeError or KeyError:
            pass


        pcs = [
            np.array([[p_cure0_u, p_cure0_dr_u, 0]]).T,
            np.array([[p_cure0_e, p_cure0_dr_e, 0]]).T,
            np.array([[p_cure0_i, 0, 0]]).T
        ]

        tps = np.zeros(3)
        for sector in range(3):
            tp0 = r_csi * pdx[sector] * y[I.Sym]
            tp1 = r_recsi * pdx[sector] * y[I.ExCS]
            tpr = r_recsi * pdx[sector] * y[I.ReCS]

            tp = tp0 + tp1 + tpr
            tps[sector] = tp.sum()

            dy[I.Sym] -= tp0
            dy[I.ExCS] -= tp1
            dy[I.ReCS] -= tpr
            p_det = np.array([1, p_dx_xpert[sector], p_dx_xpert[sector]]).reshape((3, 1))

            for sector_tx, tar_tx in zip([0, 1, 2], [I.TxPub, I.TxPriOnPub, I.TxPriOnPri]):
                txi = tp * p_det * alo[sector, sector_tx] * pcs[sector_tx]
                txf = tp * p_det * alo[sector, sector_tx] * (1 - pcs[sector_tx])
                dy[I.ReCS] += txf
                dy[tar_tx] += txi

        fn = r_csi * p_fn
        dy[I.ExCS] -= fn
        dy[I.ReCS] += fn

        fps = tps * (1 - ppv) / ppv

        da[I.A_NotiPub] += tps[0] + fps[0]
        da[I.A_NotiPri] += tps[1] + fps[1]
        return dy
