from scipy.interpolate import interp1d
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tx']


class IntvTx:
    def __init__(self, pp, uptake_time, uptake_pub, uptake_eng, rel_pltfu=0.5, p_rel=0.9):
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
        self.RateRelapse = p_rel / (1 - p_rel)

    def modify_txi(self, t, p_txi):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        if wtu <= 0 and wte <= 0:
            return p_txi

        p_txi1 = p_txi.copy()
        p_txi1[0] = 1 - (1 - p_txi1[0]) * (1 - wtu * self.RelPLTFU)
        p_txi1[1] = 1 - (1 - p_txi1[1]) * (1 - wte * self.RelPLTFU)
        return p_txi1

    def modify_rel(self, t, r_rel_teu, r_rel_tei):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)

        if wtu > 0 and r_rel_teu > self.RateRelapse:
            r_rel_teu1 = (1 - wtu) * r_rel_teu + wtu * self.RateRelapse
        else:
            r_rel_teu1 = r_rel_teu

        if wte > 0 and r_rel_tei > self.RateRelapse:
            r_rel_tei1 = (1 - wte) * r_rel_tei + wte * self.RateRelapse

        else:
            r_rel_tei1 = r_rel_tei

        return r_rel_teu1, r_rel_tei1


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
        return IntvTx(p, x, y_pub, y_eng, 0.5,  p_rel=0.1)
    elif key == 'LA-INJ':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 4 + [.10, .30, .60, .90, .90]
        y_eng = [0] * 5 + [.20, .40, .60, .60]
        return IntvTx(p, x, y_pub, y_eng, 0.25, p_rel=0.1)
    elif key == 'NoRelapse':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(p, x, y_pub, y_eng, 0, p_rel=0)
    elif key == 'NoPLTFU':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.5, 1.0] + [1.0] * 6
        return IntvTx(p, x, y_pub, y_eng, 1)
    elif key == 'HighPPM':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0, 0.5, 1.0] + [1.0] * 6
        y_eng = [0, 0.375, 0.75] + [0.75] * 6
        return IntvTx(p, x, y_pub, y_eng, 0, p_rel=0.05)

    return None

