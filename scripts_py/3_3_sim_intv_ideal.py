__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    import numpy as np
    from dy.intervention import compose_intv
    from sim.dy.obj import load_obj
    import pandas as pd
    from tqdm import tqdm
    import json

    folder = '../out/dyage'

    obj = load_obj(
        folder_input='../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        suffix='cas_cdx',
        agp=True
    )

    with open(f'{folder}/Sim_Baseline.json', 'r') as f:
        y0s = json.load(f)

    mss = list()

    for i, state in tqdm(enumerate(y0s[:100])):
        y0, pars = state['Y0'], state['Pars']
        y0 = np.array(y0)

        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)

        intvs = {
            'Baseline': compose_intv(p),
            'PerfectDx': compose_intv(p, dx='PerfectDx'),
            'PerfectTxi': compose_intv(p, dx='PerfectTxi'),
            'NoRelapse': compose_intv(p, tx='NoRelapse'),
        }

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'{folder}/Sim_IntvIdeal.csv')
