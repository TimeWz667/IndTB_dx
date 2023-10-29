import pandas as pd
from sims_pars import bayes_net_from_script
from sims_pars.fit.targets import read_targets
from sims_pars.fit.base import DataModel, Particle
from sim.inputs import load_inputs
from sim.dy import ModelBaseline #, ModelAgeGrp
import numpy as np
from scipy.stats import binom

__all__ = ['load_obj_baseline', 'load_obj_age']


class Obj(DataModel):
    def __init__(self, model, file_prior, tars, exo=None, agp=False):
        with open(file_prior, 'r') as f:
            scr = f.read()
        bn = bayes_net_from_script(scr)

        dat = self.read_data(tars, agp=agp)

        DataModel.__init__(self, dat, bn, exo=exo)

        self.Model = model
        self.Cas = self.Model.Inputs.Cascade

    @staticmethod
    def read_data(tars, agp):
        tar_inc = tars[tars.Index == 'IncR']
        tar_inc = tar_inc[tar_inc.Tag == 'All']
        tar_inc = tar_inc[tar_inc.Year >= 2015]

        tar_prev = tars[tars.Index == 'PrevUt']

        if agp:
            tar_inc_age = tars[tars.Index == 'IncR']
            tar_inc_age = tar_inc_age[tar_inc_age.Tag in ['15-24', '25-34', '35-44', '45-54', '55-64', '65+']]
            tars = pd.concat([tar_inc, tar_prev, tar_inc_age])
        else:
            tars = pd.concat([tar_inc, tar_prev])

        tars = {f'{row.Index}_{row.Tag}_{row.Year:d}': dict(row) for i, row in tars.iterrows()}
        # tars['PrPub|Treated'] = {'N': 307, 'M': 0.82}
        tars['PrevDR'] = {'N': 2404291, 'M': 0.027}

        for d in tars.values():
            d['m'] = d['M']
            d['n'] = d['N']
            d['l'], d['u'] = binom.interval(0.95, n=np.round(d['n']), p=d['M'])
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

        for t in range(2015, 2022):
            ext[f'IncR_All_{t:d}'] = ms.IncR[t]
        ext['PrevUt_All_2019'] = ms.PrevUt[2020]

        # ext['PrPub|Treated'] = ms.IncTreatedPubR[2020] / ms.IncTreatedR[2020]
        ext['PrevDR'] = ms.PrevDR[2021]

        if 'IncR_25-34' in ms:
            for k in ['15-24', '25-34', '35-44', '45-54', '55-64', '65+']:
                ext[f'IncR_{k}_2021'] = ms[f'IncR_{k}'][2021]

        return ext

    def simulate(self, pars) -> Particle:
        p = self.Cas.prepare_pars(pars)
        _, ms, _ = self.Model.simulate_to_fit(p, t_eval=np.linspace(2014, 2021, 8))
        ext = self.map_sim_data(ms)
        return Particle(pars, ext)


def load_obj_baseline(folder_input, file_prior, file_targets, year0=2000, exo=None, suffix='cas_cdx'):
    inp = load_inputs(folder_input, cs_suffix=suffix)
    inp.Demography.set_year0(year0)
    model = ModelBaseline(inp)
    targets = pd.read_csv(file_targets)
    return Obj(model, file_prior, tars=targets, exo=exo)


def load_obj_age(folder_input, file_prior, file_targets, year0=2000, exo=None, suffix='cas_cdx'):
    inp = load_inputs(folder_input, cs_suffix=suffix, agp='who')
    inp.Demography.set_year0(year0)
    model = ModelAgeGrp(inp)
    targets = pd.read_csv(file_targets)
    return Obj(model, file_prior, tars=targets, exo=exo)


if __name__ == '__main__':
    exo = {
        'drt_act': 0,
        'drt_trans': 0
    }

    obj = load_obj_baseline(
        folder_input=f'../../pars',
        file_prior='../../data/prior.txt',
        file_targets='../../data/targets.csv',
        year0=2000,
        exo=exo
    )

    # Test Objectives
    p1 = obj.sample_prior()
    p1['beta'] = 15
    tofit = obj.simulate(p1)
    print(tofit.Sims)
    print(obj.calc_distance(tofit))

    print('---------' + 'Age Group' + '_' * 10)

    # obj = load_obj_age(
    #     folder_input=f'../../pars',
    #     file_prior='../../data/prior.txt',
    #     file_targets='../../data/targets.csv',
    #     year0=2000,
    #     exo=exo
    # )
    #
    # # Test Objectives
    # p1 = obj.sample_prior()
    # p1['beta'] = 15
    # tofit = obj.simulate(p1)
    # print(tofit.Sims)
    # print(obj.calc_distance(tofit))
