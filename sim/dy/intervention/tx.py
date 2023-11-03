from scipy.interpolate import interp1d
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_intv_tx', 'IntvTx']


class IntvTx:
    def __init__(self, uptake_time, uptake_pub, uptake_eng, p_cure=0.85, p_cure_sl=0.3, rel_pltfu=0.5, p_rel=0.9):
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
        self.PCure = p_cure
        self.PCureSl = p_cure_sl
        self.RelPLTFU = rel_pltfu
        self.RateRelapse = p_rel / (1 - p_rel)

    def modify_txi(self, t, p_txi):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)
        if wtu <= 0 and wte <= 0:
            return p_txi

        p_txi1 = p_txi.copy()
        p_txi1[0] = p_txi1[0] + (1 - p_txi1[0]) * wtu * self.RelPLTFU
        p_txi1[1] = p_txi1[1] + (1 - p_txi1[1]) * wte * self.RelPLTFU
        return p_txi1

    def modify_cure(self, t, p_cure_u, p_cure_e, p_cure_dr_u, p_cure_dr_e):
        wtu, wte = self.UptakePub(t), self.UptakeEng(t)

        if wtu > 0:
            if self.PCure > p_cure_u:
                p_cure_u1 = self.PCure * wtu + p_cure_u * (1 - wtu)
            else:
                p_cure_u1 = p_cure_u
            if self.PCureSl > p_cure_dr_u:
                p_cure_dr_u1 = self.PCureSl * wtu + p_cure_dr_u * (1 - wtu)
            else:
                p_cure_dr_u1 = p_cure_dr_u
        else:
            p_cure_u1, p_cure_dr_u1 = p_cure_u, p_cure_dr_u

        if wte > 0:
            if self.PCure > p_cure_e:
                p_cure_e1 = self.PCure * wte + p_cure_e * (1 - wte)
            else:
                p_cure_e1 = p_cure_e
            if self.PCureSl > p_cure_dr_e:
                p_cure_dr_e1 = self.PCureSl * wte + p_cure_dr_e * (1 - wte)
            else:
                p_cure_dr_e1 = p_cure_dr_e
        else:
            p_cure_e1, p_cure_dr_e1 = p_cure_e, p_cure_dr_e

        return p_cure_u1, p_cure_e1, p_cure_dr_u1, p_cure_dr_e1

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


def get_intv_tx(key):
    if key == 'BPaLM':
        x = np.linspace(2024, 2040, 9)
        y_pub = [20 / 100, 80 / 100] + [1] * 7
        y_eng = [10 / 100, 60 / 100] + [1] * 7
        return IntvTx(x, y_pub, y_eng, 0, 0.85, rel_pltfu=0, p_rel=0.3)
    elif key == 'PAN-TB':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 3 + [.10, .30, .90, .90, .90, .90]
        y_eng = [0] * 4 + [.20, .60, .90, .90, .90]
        return IntvTx(x, y_pub, y_eng, 0.92, 0, rel_pltfu=0.5, p_rel=0.12)
    elif key == 'LA-INJ':
        x = np.linspace(2024, 2040, 9)
        y_pub = [0] * 4 + [.10, .30, .60, .90, .90]
        y_eng = [0] * 5 + [.20, .40, .60, .60]
        return IntvTx(x, y_pub, y_eng, 0.98, 0, rel_pltfu=1, p_rel=0.05)
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
        return IntvTx(x, y_pub, y_eng, p_cure=0.9, p_cure_sl=0, rel_pltfu=0.5, p_rel=s)
    return None

