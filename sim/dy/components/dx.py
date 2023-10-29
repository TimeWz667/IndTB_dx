from sim.dy.components.proc import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Dx']


class Dx(Process):
    def __init__(self, keys):
        Process.__init__(self, 'dx', keys)

    def calculate_calc(self, t, y, pars, calc: dict, **kwargs):
        I = self.Keys

        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']

        k = np.exp(pars['rt_cs'] * max(t - 2019, pars['t0_decline'] - 2019))
        r_csi *= k
        r_recsi *= k

        p_ent, p_itt, p_pdx, p_txi = pars['p_ent'], pars['p_itt'], pars['p_dx'], pars['p_txi']
        p_dx_xpert = pars['p_dx_xpert']

        try:
            intv_dx = kwargs['intv'].Dx
            p_ent, p_pdx, p_txi = intv_dx.modify_dx(t, p_ent, p_pdx, p_txi)
        except AttributeError or KeyError:
            pass

        try:
            intv_ppm = kwargs['intv'].PPM
            p_ent = intv_ppm.modify_ent(t, p_ent)
        except AttributeError or KeyError:
            pass

        try:
            intv_tx = kwargs['intv'].Tx
            p_txi = intv_tx.modify_txi(t, p_txi)
        except AttributeError or KeyError:
            pass

        pdx = p_ent * p_itt * p_pdx * p_txi
        ppv, alo = pars['ppv'], pars['tx_alo']
        p_cure0, p_cure0_dr = pars['p_cure'], pars['p_cure_dr']
        p_cure0_u, p_cure0_e, p_cure0_i = p_cure0, p_cure0, p_cure0
        p_cure0_dr_u, p_cure0_dr_e, p_cure0_dr_i = p_cure0_dr, p_cure0_dr, p_cure0_dr

        try:
            intv_tx = kwargs['intv'].Tx
            p_cure0_u, p_cure0_e, p_cure0_dr_u, p_cure0_dr_e = intv_tx.modify_cure(t, p_cure0_u, p_cure0_e,
                                                                                   p_cure0_dr_u, p_cure0_dr_e)
        except AttributeError or KeyError:
            pass

        pcs = [
            np.array([p_cure0_u, p_cure0_dr_u, p_cure0_dr_u]).reshape((1, 1, 3)),
            np.array([p_cure0_e, p_cure0_dr_e, p_cure0_dr_e]).reshape((1, 1, 3)),
            np.array([p_cure0_i, p_cure0_dr_i, p_cure0_dr_i]).reshape((1, 1, 3))
        ]

        calc['fp'] = fps = np.zeros_like(ppv)

        calc['dx_pub'] = dx = {
            'tp0': r_csi * pdx[0] * y[I.Sym],
            'fn0': r_csi * (p_ent[0] - pdx[0]) * y[I.Sym],
            'tp1': r_recsi * pdx[0] * y[I.ExCS],
            'tpr': r_recsi * pdx[0] * y[I.ReCS]
        }
        dx['tp'] = tp = dx['tp0'] + dx['tp1'] + dx['tpr']

        p_cure = pcs[0] * np.array([1, p_dx_xpert[0], p_dx_xpert[0]]).reshape((1, 1, 3))
        dx['txi'] = tp.reshape((-1, *tp.shape)) * alo[0].reshape((-1, 1, 1)) * p_cure
        dx['txf'] = tp.reshape((-1, *tp.shape)) * alo[0].reshape((-1, 1, 1)) * (1 - p_cure)

        fps[0] = tp.sum() * (1 - ppv[0]) / ppv[0]

        calc['dx_eng'] = dx = {
            'tp0': r_csi * pdx[1] * y[I.Sym],
            'fn0': r_csi * (p_ent[1] - pdx[1]) * y[I.Sym],
            'tp1': r_recsi * pdx[1] * y[I.ExCS],
            'tpr': r_recsi * pdx[1] * y[I.ReCS]
        }
        dx['tp'] = tp = dx['tp0'] + dx['tp1'] + dx['tpr']

        p_cure = pcs[0] * np.array([1, p_dx_xpert[1], p_dx_xpert[1]]).reshape((1, 1, 3))
        dx['txi'] = tp.reshape((-1, *tp.shape)) * alo[1].reshape((-1, 1, 1)) * p_cure
        dx['txf'] = tp.reshape((-1, *tp.shape)) * alo[1].reshape((-1, 1, 1)) * (1 - p_cure)

        fps[1] = tp.sum() * (1 - ppv[1]) / ppv[1]

        calc['dx_pri'] = dx = {
            'tp0': r_csi * pdx[2] * y[I.Sym],
            'fn0': r_csi * (p_ent[2] - pdx[2]) * y[I.Sym],
            'tp1': r_recsi * pdx[2] * y[I.ExCS],
            'tpr': r_recsi * pdx[2] * y[I.ReCS]
        }
        dx['tp'] = tp = dx['tp0'] + dx['tp1'] + dx['tpr']

        p_cure = pcs[0] * np.array([1, p_dx_xpert[2], p_dx_xpert[2]]).reshape((1, 1, 3))
        dx['txi'] = tp.reshape((-1, *tp.shape)) * alo[2].reshape((-1, 1, 1)) * p_cure
        dx['txf'] = tp.reshape((-1, *tp.shape)) * alo[2].reshape((-1, 1, 1)) * (1 - p_cure)

        fps[2] = tp.sum() * (1 - ppv[2]) / ppv[2]

        return calc

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        for dx in [calc['dx_pub'], calc['dx_eng'], calc['dx_pri']]:
            tp0, fn0, tp1, tpr = dx['tp0'], dx['fn0'], dx['tp1'], dx['tpr']
            txi, txf = dx['txi'], dx['txf']

            dy[I.Sym] += - tp0 - fn0
            dy[I.ExCS] += fn0 - tp1
            dy[I.ReCS] += txf.sum(0) - tpr
            dy[I.TxPub] += txi[0]
            dy[I.TxPriOnPub] += txi[1]
            dy[I.TxPriOnPri] += txi[2]

        fp = calc['fp']

        da[I.A_NotiPub] += (calc['dx_pub']['tp0'] + calc['dx_pub']['tp1']).sum() + fp[0].sum()
        da[I.A_NotiPri] += (calc['dx_eng']['tp0'] + calc['dx_eng']['tp1']).sum() + fp[1].sum()

        return dy, da

    def measure(self, mea, t, y, pars, calc, **kwargs):
        n = mea['N']

        mea.update({
            'Tp_pub': calc['tp0'][0] + calc['tp1'][0],
            'Tp_eng': calc['tp0'][1] + calc['tp1'][1],
            'Tp_pri': calc['tp0'][2] + calc['tp1'][2],
            'Fp_pub': calc['fp'][0],
            'Fp_eng': calc['fp'][1],
            'Fp_pri': calc['fp'][2]
        })
        mea['CNR_Pub'] = (mea['Tp_pub'] + mea['Fp_pub']) / n
        mea['CNR_Pri'] = (mea['Tp_eng'] + mea['Fp_eng']) / n
        mea['CNR'] = mea['CNR_Pub'] + mea['CNR_Pri']


if __name__ == '__main__':
    from sim.inputs import load_inputs

    inp = load_inputs('../../../pars', cs_suffix='bac_cdx_sector_2022')

    print(inp)

