from scipy.interpolate import interp1d
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tx']


class IntvTx:
    def __init__(self, pp, uptake_time, uptake_pub, uptake_eng, rel_pltfu=0.5, txs=0.95, rel_rel=1.0):
        self.UptakePub = interp1d(
            x=np.concatenate([[0, 2023], uptake_time]),
            y=np.concatenate([[0, 0], uptake_pub]),
            fill_value="extrapolate"
        )
        self.UptakeEng = interp1d(
            x=np.concatenate([[0, 2023], uptake_time]),
            y=np.concatenate([[0, 0], uptake_eng]),
            fill_value="extrapolate"
        )
        self.RelPLTFU = rel_pltfu
        self.PrTxs = txs
        self.RelRelapse = rel_rel

    def modify_txi(self, t, p_txi):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        if wtu <= 0 and wte <= 0:
            return p_txi

        p_txi1 = p_txi.copy()
        p_txi1[0] = 1 - (1 - p_txi1[0]) * (1 - wtu * self.RelPLTFU)
        p_txi1[1] = 1 - (1 - p_txi1[1]) * (1 - wte * self.RelPLTFU)
        return p_txi1

    def modify_txo(self, t, r_txs, r_txl):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        if wtu <= 0 and wte <= 0:
            return r_txs, r_txl
        p0u, p0e = (r_txs / (r_txl + r_txs))[:2]
        p1u, p1e = p0u + wtu * (self.PrTxs - p0u), p0e + wte * (self.PrTxs - p0e)
        p1u, p1e = max(p1u, p0u), max(p1e, p0e)

        r_txs1, r_txl1 = r_txs.copy(), r_txl.copy()
        r_txs1[0], r_txl1[0] = p1u * (r_txs1[0] + r_txl1[0]), (1 - p1u) * (r_txs1[0] + r_txl1[0])
        r_txs1[1], r_txl1[1] = p1e * (r_txs1[1] + r_txl1[1]), (1 - p1e) * (r_txs1[1] + r_txl1[1])
        return r_txs1, r_txl1

    def modify_relu(self, t, r_rel_tcu, r_rel_tdu, r_rel_stu):
        wt = self.UptakePub(t)
        if wt <= 0 or self.RelRelapse >= 1:
            return r_rel_tcu, r_rel_tdu, r_rel_stu

        rel = (1 - wt) + wt * self.RelRelapse

        r_rel_tcu, r_rel_tdu = r_rel_tcu * rel, r_rel_tdu * rel
        r_rel_tcu, r_rel_tdu = max(r_rel_tcu, r_rel_stu), max(r_rel_tdu, r_rel_stu)
        return r_rel_tcu, r_rel_tdu, r_rel_stu

    def modify_reli(self, t, r_rel_tci, r_rel_tdi, r_rel_sti):
        wt = self.UptakeEng(t)
        if wt <= 0 or self.RelRelapse >= 1:
            return r_rel_tci, r_rel_tdi, r_rel_sti

        rel = (1 - wt) + wt * self.RelRelapse

        r_rel_tci, r_rel_tdi = r_rel_tci * rel, r_rel_tdi * rel
        r_rel_tci, r_rel_tdi = max(r_rel_tci, r_rel_sti), max(r_rel_tdi, r_rel_sti)
        return r_rel_tci, r_rel_tdi, r_rel_sti


def get_intv_tx(p, key):
    if key == 'BPaLM':
        x = np.linspace(2024, 2040, 9)
        y_pub = [20 / 100, 80 / 100] + [1] * 7
        y_eng = [10 / 100, 60 / 100] + [1] * 7

        return None # todo
    elif key == 'PAN-TB':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 3 + [.10, .30, .90, .90, .90, .90]
        y_eng = [0] * 4 + [.20, .60, .90, .90, .90]
        return IntvTx(p, x, y_pub, y_eng, 0.5, 0.9, rel_rel=0.1)
    elif key == 'LA-INJ':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 4 + [.10, .30, .60, .90, .90]
        y_eng = [0] * 5 + [.20, .40, .60, .60]
        return IntvTx(p, x, y_pub, y_eng, 0.25, 0.9, rel_rel=0.1)
    elif key == 'NoRelapse':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(p, x, y_pub, y_eng, 0, 0.9, 0)
    elif key == 'NoPLTFU':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(p, x, y_pub, y_eng, 1, 0.9, 0)

    return None

