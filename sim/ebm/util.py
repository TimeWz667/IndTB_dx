from abc import ABCMeta, abstractmethod
from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd

__author__ = 'Chu-Chang Ku'
__all__ = ['AbsModelODE']


def diff_mea(mss0, mss1, dt):
    keys = list(mss0.columns)
    mss1.index = mss0.index

    mss = {k: (mss0[k] + mss1[k]) / 2 for k in keys}
    n = mss['N']

    for k in keys:
        if k.startswith('Cum'):
            mss[k.replace('Cum', '') + 'R'] = (mss1[k] - mss0[k]) / n * dt

    mss = pd.DataFrame.from_dict(mss)
    return mss


class AbsModelODE(metaclass=ABCMeta):
    def __init__(self, n_dim, inputs, t0, t1, dt, t_warmup=0, dfe=None):
        self.NDim = n_dim
        self.T0 = t0
        self.T1 = t1
        self.DT = dt
        self.TWarmup = t_warmup
        self.DFE = dfe
        self.Inputs = inputs

    @abstractmethod
    def reform_parameters(self, p):
        pass

    @abstractmethod
    def get_y0(self, pars) -> np.ndarray:
        pass

    @abstractmethod
    def __call__(self, t, y, pars, intv=None):
        pass

    @abstractmethod
    def measure(self, t, y, pars, intv=None) -> dict:
        pass

    def _simulate(self, t0, t1, y0, pars, intv=None, dense_output=False):
        y0 = np.array(y0).reshape(-1)
        pars = self.reform_parameters(pars)

        if intv is None:
            ys = solve_ivp(self, [t0, t1], y0, args=(pars, ), events=self.DFE, dense_output=dense_output)
        else:
            ys = solve_ivp(self, [t0, t1], y0, args=(pars, intv, ), events=self.DFE, dense_output=dense_output)

        if self.DFE is not None:
            if len(ys.t_events[0]) > 0 or not ys.success:
                return None, None, {'succ': False, 'res': 'DFE reached'}

        return ys, {'succ': True, 'TimeSpan': [t0, t1], 'y': ys.y.T[-1]}

    def simulate_to_baseline(self, pars):
        y0 = self.get_y0(pars)
        ys, msg = self._simulate(self.T0 - self.TWarmup, self.T0, y0=y0, pars=pars,
                                 dense_output=False)
        if not msg['succ']:
            return None, None, msg

        return ys, None, msg

    def simulate_onward(self, y0, pars, intv=None):
        t0, t1 = self.T0, self.T1 + self.DT

        pars = self.reform_parameters(pars)

        ys, msg = self._simulate(t0, t1, y0=y0, pars=pars, intv=intv, dense_output=True)
        t_eval = np.linspace(t0, t1, int((t1 - t0) / self.DT + 1))

        if msg['succ']:
            t_eval0, t_eval1 = t_eval, t_eval + self.DT
            t_all = list(set(t_eval0).union(set(t_eval1)))
            t_all.sort()

            mss = {t: self.measure(t, ys.sol(t), pars, intv=intv) for t in t_all}
            mss0 = pd.DataFrame([mss[t] for t in t_eval0]).set_index('Year')
            mss1 = pd.DataFrame([mss[t] for t in t_eval1]).set_index('Year')
            mss = diff_mea(mss0, mss1, self.DT)

            return ys, mss, msg
        else:
            return None, None, msg

    def simulate(self, pars, intv=None):
        _, _, msg = self.simulate_to_baseline(pars)
        assert msg['succ'] is True
        y0 = msg['y']
        ys, mss, msg = self.simulate_onward(y0, pars=pars, intv=intv)
        assert msg['succ'] is True
        return ys, mss, msg

    def simulate_to_fit(self, pars, t_eval, **kwargs):
        y0 = self.get_y0(pars)
        t0 = min(self.T0, min(t_eval))

        ys, msg = self._simulate(t0 - self.TWarmup, t0, y0=y0, pars=pars,
                                 dense_output=False, **kwargs)
        if not msg['succ']:
            return None, None, msg

        y0 = msg['y']
        t1 = max(t_eval) + self.DT
        ys, msg = self._simulate(t0, t1, y0=y0, pars=pars, dense_output=True, **kwargs)

        if msg['succ']:
            p = self.reform_parameters(pars)

            t_eval0, t_eval1 = t_eval, t_eval + self.DT
            t_all = list(set(t_eval0).union(set(t_eval1)))
            t_all.sort()

            mss = {t: self.measure(t, ys.sol(t), p, **kwargs) for t in t_all}
            mss0 = pd.DataFrame([mss[t] for t in t_eval0]).set_index('Year')
            mss1 = pd.DataFrame([mss[t] for t in t_eval1]).set_index('Year')
            mss = diff_mea(mss0, mss1, self.DT)

            return ys, mss, msg
        else:
            return None, None, msg
