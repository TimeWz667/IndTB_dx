from sims_pars import bayes_net_from_script
from sims_pars.fit.targets import read_targets
from sims_pars.fit.base import DataModel, Particle
from sim.inputs import load_inputs
from sim.ebm.dy import ModelBaseline
import numpy as np
from scipy.stats import binom

__all__ = ['load_obj_baseline']


class Obj(DataModel):
    def __init__(self, model, file_prior, tars, exo=None):
        with open(file_prior, 'r') as f:
            scr = f.read()
        bn = bayes_net_from_script(scr)

        dat = self.read_data(tars)

        DataModel.__init__(self, dat, bn, exo=exo)

        self.Model = model
        self.Cas = self.Model.Inputs.Cascade

    @staticmethod
    def read_data(tars):
        for d in tars.values():
            d['l'], d['u'] = binom.interval(0.95, n=np.round(d['n']), p=d['m'])
            d['l'] /= d['n']
            d['u'] /= d['n']

        dat = read_targets(tars)
        for d in dat.values():
            try:
                d.Range = np.abs(d.Range)
            except AttributeError:
                pass
        return dat

    @staticmethod
    def map_sim_data(ms):
        ext = dict()

        for t in range(2018, 2022):
            ext[f'Det_Pub_{t:d}'] = ms.NotiPubR[t]
            ext[f'Det_Eng_{t:d}'] = ms.NotiPriR[t]
            # ext[f'CNR_All_{t:d}'] = ms.CNR[t]
        ext['Prev_Ut_2020'] = ms.PrevUt[2020]
        ext['Prev_Asym_2020'] = ms.PrevA[2020]
        ext['Prev_Sym_2020'] = ms.PrevS[2020]
        ext['Prev_ExCS_2020'] = ms.PrevC[2020]
        # ext['PrA_All_2019'] = ms.PrA[2019]
        # ext['PrS_All_2019'] = ms.PrS[2019]
        # ext['PrC_All_2019'] = ms.PrC[2019]

        return ext

    def simulate(self, pars) -> Particle:
        p = self.Cas.prepare_pars(pars)
        _, ms, _ = self.Model.simulate_to_fit(p, t_eval=np.linspace(2014, 2021, 8))
        ext = self.map_sim_data(ms)
        return Particle(pars, ext)


def load_obj_baseline(folder_input, file_prior, year0=2000, exo=None):
    inp = load_inputs(folder_input)
    inp.Demography.set_year0(year0)
    model = ModelBaseline(inp)
    return Obj(model, file_prior, tars=inp.Targets, exo=exo)
