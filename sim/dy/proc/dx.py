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

        # Diagnosis
        k = np.exp(pars['rt_cs'] * max(t - pars['t0_decline'], 0))
        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']
        r_csi *= k
        r_recsi *= k
        ppv, alo = pars['ppv'], pars['tx_alo']

        p_ent, p_itt, p_pdx, p_txi = pars['p_ent'], pars['p_itt'], pars['p_dx'], pars['p_txi']
        p_txi_new = np.zeros_like(p_txi)
        p_dx_xpert = pars['p_dx_xpert']

        try:
            p_itt, p_ent, p_pdx, p_txi = intv.Dx.modify_dx(t, p_itt, p_ent, p_pdx, p_txi)
        except (AttributeError, TypeError):
            pass

        try:
            p_ent = intv.PPM.modify_ent(t, p_ent)
        except (AttributeError, TypeError):
            pass

        try:
            p_txi_new = intv.Tx.get_p_txi_new(t, p_txi)
        except (AttributeError, TypeError):
            p_txi_new = p_txi.copy()

        try:
            p_up_soc, p_up_new = intv.Tx.modify_uptake_new(t)
        except (AttributeError, TypeError):
            p_up_soc, p_up_new = np.ones(3), np.zeros(3)

        pdx = p_ent * p_itt * p_pdx * p_up_soc * p_txi
        pdx_new = p_ent * p_itt * p_pdx * p_up_new * p_txi_new
        p_fn = 1 - pdx.sum() - pdx_new.sum()
        p_cure0, p_cure0_dr = pars['p_cure'], pars['p_cure_dr']
        p_cure0_u, p_cure0_e, p_cure0_i = p_cure0, p_cure0, p_cure0
        p_cure0_dr_u, p_cure0_dr_e, p_cure0_dr_i = p_cure0_dr, p_cure0_dr, p_cure0_dr

        try:
            p_cure0_dr_u, p_cure0_dr_e = intv.Txs.modify_cure_sl(t, p_cure0_dr_u, p_cure0_dr_e)
        except (AttributeError, TypeError):
            pass

        pcs = [
            np.array([[p_cure0_u, p_cure0_dr_u, 0]]).T,
            np.array([[p_cure0_e, p_cure0_dr_e, 0]]).T,
            np.array([[p_cure0_i, 0, 0]]).T
        ]

        try:
            p_cure0_u, p_cure0_e = intv.Tx.get_p_cure_new(p_cure0_u, p_cure0_e)
        except AttributeError or KeyError:
            pass

        pcs_new = [
            np.array([[p_cure0_u, p_cure0_dr_u, 0]]).T,
            np.array([[p_cure0_e, p_cure0_dr_e, 0]]).T,
            np.array([[p_cure0_i, 0, 0]]).T
        ]

        tps = np.zeros(3)
        for sector in range(3):
            p_det = np.array([1, p_dx_xpert[sector], p_dx_xpert[sector]]).reshape((3, 1))

            tp0 = r_csi * pdx[sector] * y[I.Sym]
            tp1 = r_recsi * pdx[sector] * y[I.ExCS]
            tpr = r_recsi * pdx[sector] * y[I.ReCS]

            tp = tp0 + tp1 + tpr
            tps[sector] += tp.sum()

            dy[I.Sym] -= tp0
            dy[I.ExCS] -= tp1
            dy[I.ReCS] -= tpr

            for sector_tx, tar_tx in zip([0, 1, 2], [I.TxPub, I.TxPriOnPub, I.TxPriOnPri]):
                txi = tp * p_det * alo[sector, sector_tx] * pcs[sector_tx]
                txf = tp * p_det * alo[sector, sector_tx] * (1 - pcs[sector_tx])
                dy[I.ReCS] += txf
                dy[tar_tx] += txi

            tp0 = r_csi * pdx_new[sector] * y[I.Sym]
            tp1 = r_recsi * pdx_new[sector] * y[I.ExCS]
            tpr = r_recsi * pdx_new[sector] * y[I.ReCS]

            tp = tp0 + tp1 + tpr
            tps[sector] += tp.sum()

            dy[I.Sym] -= tp0
            dy[I.ExCS] -= tp1
            dy[I.ReCS] -= tpr

            for sector_tx, tar_tx in zip([0, 1, 2], [I.TxNewPub, I.TxNewPriOnPub, I.TxNewPriOnPri]):
                txi = tp * p_det * alo[sector, sector_tx] * pcs_new[sector_tx]
                txf = tp * p_det * alo[sector, sector_tx] * (1 - pcs_new[sector_tx])
                dy[I.ReCS] += txf
                dy[tar_tx] += txi

        fn = r_csi * p_fn
        dy[I.Sym] -= fn
        dy[I.ExCS] += fn

        # Tx outcome
        r_txo = 1 / pars['tx_dur']
        txo = y[[I.TxPub, I.TxPriOnPub, I.TxPriOnPri]] * r_txo.reshape((-1, 1, 1))

        calc['Txo_soc'] = txo.sum((1, 2))

        try:
            r_txo_new = intv.Tx.get_r_txo(r_txo)
        except AttributeError or KeyError:
            r_txo_new = 1 / pars['tx_dur']

        txo_new = y[[I.TxNewPub, I.TxNewPriOnPub, I.TxNewPriOnPri]] * r_txo_new.reshape((-1, 1, 1))

        calc['Txo_new'] = txo_new.sum((1, 2))

        dy[I.TxPub] -= txo[0]
        dy[I.TxPriOnPub] -= txo[1]
        dy[I.TxPriOnPri] -= txo[2]

        dy[I.TxNewPub] -= txo_new[0]
        dy[I.TxNewPriOnPub] -= txo_new[1]
        dy[I.TxNewPriOnPri] -= txo_new[2]

        dy[I.RHighPub] += txo[0] + txo_new[0]
        dy[I.RHighPri] += txo[1] + txo[2] + txo_new[1] + txo_new[2]

        # ACF
        try:
            r_acf, r_acf_new = intv.ACF.modify_acf(t, 0)
        except (TypeError, AttributeError):
            r_acf = r_acf_new = 0

        acf_a = r_acf * y[I.Asym]
        acf_s = r_acf * y[I.Sym]
        acf_c = r_acf * y[I.ExCS]
        acf_r = r_acf * y[I.ReCS]
        calc['acf'] = acf = (acf_a + acf_s + acf_c + acf_r).sum()

        dy[I.Asym] -= acf_a
        dy[I.Sym] -= acf_s
        dy[I.ExCS] -= acf_c
        dy[I.ReCS] -= acf_r
        dy[I.TxPub] += acf_a + acf_s + acf_c + acf_r
        da[I.A_ACF] += acf

        acf_a = r_acf_new * y[I.Asym]
        acf_s = r_acf_new * y[I.Sym]
        acf_c = r_acf_new * y[I.ExCS]
        acf_r = r_acf_new * y[I.ReCS]
        acf = (acf_a + acf_s + acf_c + acf_r).sum()
        calc['acf'] += acf

        dy[I.Asym] -= acf_a
        dy[I.Sym] -= acf_s
        dy[I.ExCS] -= acf_c
        dy[I.ReCS] -= acf_r
        dy[I.TxNewPub] += acf_a + acf_s + acf_c + acf_r
        da[I.A_ACF] += acf

        fps = tps * (1 - ppv) / ppv

        da[I.A_NotiPub] += tps[0] + fps[0]
        da[I.A_NotiPri] += tps[1] + fps[1]

        return dy
