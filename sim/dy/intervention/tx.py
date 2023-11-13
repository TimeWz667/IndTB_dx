from scipy.interpolate import interp1d
import numpy as np
from sim.dy.util import scale_up

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tx', 'IntvTx', 'TxBPaLM']


class TxBPaLM:
    def __init__(self):
        x = np.linspace(2024, 2040, 9)
        y_pub = [20 / 100, 80 / 100] + [1] * 7
        y_eng = [10 / 100, 60 / 100] + [1] * 7

        self.UptakePub = interp1d(
            x=np.concatenate([[0, 2023], x]), y=np.concatenate([[0, 0], y_pub]), fill_value="extrapolate")
        self.UptakeEng = interp1d(
            x=np.concatenate([[0, 2023], x]), y=np.concatenate([[0, 0], y_eng]), fill_value="extrapolate")
        self.PCureSl = 0.85

    def modify_cure_sl(self, t, p_cure_dr_u, p_cure_dr_e):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        p_cure_dr_u1 = scale_up(p_cure_dr_u, self.PCureSl, wtu)
        p_cure_dr_e1 = scale_up(p_cure_dr_e, self.PCureSl, wte)
        return p_cure_dr_u1, p_cure_dr_e1


class IntvTx:
    def __init__(self, uptake_time, uptake_pub, uptake_eng, p_cure=0.85, rel_pltfu=0.5, p_rel=0.9, dur=0.5):
        self.UptakePub = interp1d(
            x=np.concatenate([[0, 2023], uptake_time]), y=np.concatenate([[0, 0], uptake_pub]), fill_value="extrapolate")
        self.UptakeEng = interp1d(
            x=np.concatenate([[0, 2023], uptake_time]), y=np.concatenate([[0, 0], uptake_eng]), fill_value="extrapolate")
        self.PCure = p_cure
        self.RelPLTFU = rel_pltfu
        self.RateRelapse = p_rel / (1 - p_rel)
        self.DurTx = dur

    def modify_uptake_new(self, t):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        p_up = np.array([wtu, wte, 0])
        return 1 - p_up, p_up

    def modify_cure(self, t, p_cure_u, p_cure_e):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        p_cure_u1 = scale_up(p_cure_u, self.PCure, wtu)
        p_cure_e1 = scale_up(p_cure_e, self.PCure, wte)
        return p_cure_u1, p_cure_e1

    def get_p_cure_new(self, p_cure_u, p_cure_e):
        return max(p_cure_u, self.PCure), max(p_cure_e, self.PCure)

    def get_r_rel_new(self, r_rel_u, r_rel_i):
        return min(r_rel_u, self.RateRelapse), min(r_rel_i, self.RateRelapse)

    def get_p_txi_new(self, p_txi):
        return 1 - (1 - p_txi) * self.RelPLTFU

    def get_r_txo(self, r_txo):
        return np.ones_like(r_txo) / self.DurTx


def get_intv_tx(key):
    if key == 'PAN-TB':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 3 + [.10, .30, .90, .90, .90, .90]
        y_eng = [0] * 4 + [.20, .60, .90, .90, .90]
        return IntvTx(x, y_pub, y_eng, 0.92, rel_pltfu=0.5, p_rel=0.12, dur=2 / 12)
    elif key == 'LA-INJ':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 4 + [.10, .30, .60, .90, .90]
        y_eng = [0] * 5 + [.20, .40, .60, .60]
        return IntvTx(x, y_pub, y_eng, 0.98, rel_pltfu=0.25, p_rel=0.05, dur=1 / 12)
    elif key == 'Null_PAN-TB':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 9
        y_eng = [0] * 9
        return IntvTx(x, y_pub, y_eng, 0.92, rel_pltfu=0.5, p_rel=0.12, dur=2 / 12)
    elif key == 'Null_LA-INJ':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 9
        y_eng = [0] * 9
        return IntvTx(x, y_pub, y_eng, 0.98, rel_pltfu=0.25, p_rel=0.05, dur=1 / 12)
    elif key == 'NoRelapse':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(x, y_pub, y_eng, 0, p_rel=0)
    elif key == 'NoPLTFU':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(x, y_pub, y_eng, 1)
    elif key == 'HighPPM':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.375, 0.75] + [0.75] * 6
        return IntvTx(x, y_pub, y_eng, 0, p_rel=0.05)

    elif key.startswith('RelRed_'):
        s = float(key.replace("RelRed_", ""))
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 3 + [.10, .30, .90, .90, .90, .90]
        y_eng = [0] * 4 + [.20, .60, .90, .90, .90]
        return IntvTx(x, y_pub, y_eng, p_cure=0.9, rel_pltfu=0.5, p_rel=s)
    return None

