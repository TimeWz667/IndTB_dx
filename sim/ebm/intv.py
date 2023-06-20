from pydantic import BaseModel
from pydantic.types import confloat
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class MassACF(BaseModel):
    Coverage: confloat(ge=0, le=1) = 0
    Sens_A: confloat(ge=0, le=1) = 0.98
    Sens_S: confloat(ge=0, le=1) = 0.98
    Sens_L: confloat(ge=0, le=1) = 0.4
    Spec: confloat(ge=0, le=1) = 0.6
    Target: str = 'All'  # All, 30%, 20%, 10%


class Intervention(BaseModel):
    MassACF: MassACF = MassACF()

    T0_ScaleUp: float = 2023
    T1_ScaleUp: float = 2025

    def get_wt(self, t):
        t = max(min(t, self.T1_ScaleUp), self.T0_ScaleUp)
        return (t - self.T0_ScaleUp) / (self.T1_ScaleUp - self.T0_ScaleUp)

    def modify_acf(self, t, n_dim=1):
        if t >= self.T0_ScaleUp and self.MassACF.Coverage > 0:
            acf = self.MassACF
            r_acf0 = acf.Coverage * self.get_wt(t)

            if n_dim > 1:
                eli = np.ones(4)
                if acf.Target == '10%':
                    eli[:3] = 0
                    r_acf0 = r_acf0 * eli / 0.1
                elif acf.Target == '20%':
                    eli[:2] = 0
                    r_acf0 = r_acf0 * eli / 0.2
                elif acf.Target == '30%':
                    eli[:1] = 0
                    r_acf0 = r_acf0 * eli / 0.3

            return r_acf0 * acf.Sens_A, r_acf0 * acf.Sens_S, r_acf0 * acf.Sens_L, r_acf0 * (1 - acf.Spec)
        else:
            return np.zeros(n_dim), np.zeros(n_dim), np.zeros(n_dim), np.zeros(n_dim)
