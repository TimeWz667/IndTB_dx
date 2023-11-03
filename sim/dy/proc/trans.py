import numpy as np
from sim.dy.proc.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Transmission']


class Transmission(Process):
    def __init__(self, keys):
        Process.__init__(self, 'trans', keys)

    def find_dya(self, t, y, da, pars, calc, intv=None):
        I = self.Keys

        n = y.sum()

        sus = pars['sus']
        if intv is not None:
            try:
                sus = intv.Vac.modify_sus(t, sus)
            except AttributeError or KeyError:
                pass

        dt = max(t - pars['t0_decline'], 0)
        k = np.exp(- pars['drt_trans'] * dt)

        foi_ds = (pars['trans_ds'] * y).sum() * k * pars['beta'] / n
        foi_dr = (pars['trans_dr'] * y).sum() * k * pars['beta'] * pars['rr_beta_dr'] / n
        foi_fr = (pars['trans_fr'] * y).sum() * k * pars['beta'] * pars['rr_beta_fr'] / n

        infection_ds = sus * foi_ds * y
        infection_dr = sus * foi_dr * y
        infection_fr = sus * foi_fr * y

        calc['infection_ds'] = infection_ds.sum()
        calc['infection_dr'] = infection_dr.sum()
        calc['infection_fr'] = infection_fr.sum()

        dy = - infection_ds - infection_dr - infection_fr
        dy[I.FLat, I.DS] += infection_ds.sum((0, 1))
        dy[I.FLat, I.DR] += infection_dr.sum((0, 1))
        dy[I.FLat, I.FR] += infection_fr.sum((0, 1))

        return dy
