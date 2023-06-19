from sim.components.proc import Process
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

        p_ent, p_txi = pars['p_ent'], pars['p_txi']
        pdx0 = p_ent * pars['p_dx0'] * p_txi
        pdx1 = p_ent * pars['p_dx1'] * p_txi

        calc['tp0'] = tp0 = r_csi * pdx0.reshape((-1, 1)) * y[I.Sym].reshape((1, -1))
        calc['fn0'] = r_csi * (p_ent - pdx0).reshape((-1, 1)) * y[I.Sym].reshape((1, -1))
        calc['tp1'] = tp1 = r_recsi * pdx1.reshape((-1, 1)) * y[I.ExCS].reshape((1, -1))

        tp, alo = (tp0 + tp1).reshape((3, 1, -1)), pars['tx_alo']

        calc['txi'] = (tp * alo.reshape((3, -1, 1))).sum(0)

        ppv = pars['ppv'].reshape((-1, 1))
        calc['fp'] = (tp0 + tp1) * (1 - ppv) / ppv

        if 'intv' in kwargs and kwargs['intv'] is not None:
            r_acf_a, r_acf_s, r_acf_l, r_acf_ot = kwargs['intv'].modify_acf(t, n_dim=y.shape[1])
        else:
            r_acf_a, r_acf_s, r_acf_l, r_acf_ot = 0, 0, 0, 0

        calc['acf_a'] = r_acf_a * y[I.Asym]
        calc['acf_s'] = r_acf_s * y[I.Sym]
        calc['acf_c'] = r_acf_s * y[I.ExCS]
        calc['acf_l'] = r_acf_l * y[I.FLat]
        calc['acf_td'] = r_acf_l * y[I.RHigh]

        calc['acf_ot'] = r_acf_ot * (y[I.U] + y[I.SLat] + y[I.RLow] + y[I.RSt])

        return calc

    def compose_dya(self, ya, calc: dict):
        I = self.Keys

        y, aux = ya
        dy, da = np.zeros_like(y), np.zeros_like(aux)

        tp0, fn0, tp1 = calc['tp0'], calc['fn0'], calc['tp1']

        txi = calc['txi']

        dy[I.Sym] += - tp0.sum(0) - fn0.sum(0)
        dy[I.ExCS] += fn0.sum(0) - tp1.sum(0)
        dy[I.TxPub] += txi[0]
        dy[I.TxPri] += txi[1]

        acf_a, acf_s, acf_c = calc['acf_a'], calc['acf_s'], calc['acf_c']
        dy[I.Asym] += - acf_a
        dy[I.Sym] += - acf_s
        dy[I.ExCS] += - acf_c
        dy[I.TxPub] += acf_a + acf_s + acf_c

        acf_l, acf_td = calc['acf_l'], calc['acf_td']
        dy[I.FLat] -= acf_l
        dy[I.SLat] += acf_l
        dy[I.RHigh] -= acf_td
        dy[I.RLow] += acf_td

        fp = calc['fp']

        da[I.A_NotiPub] += (tp0[0] + tp1[0] + fp[0]).sum()
        da[I.A_NotiPri] += (tp0[1] + tp1[1] + fp[1]).sum()

        acf_ot = calc['acf_ot']
        da[I.A_YieldATB] += acf_a + acf_s + acf_c
        da[I.A_YieldLTBI] += acf_l + acf_td
        da[I.A_YieldOT] += acf_ot

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
