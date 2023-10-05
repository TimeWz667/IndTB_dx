
__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    import numpy as np
    from sim.ebm.intervention import compose_intv
    from sim.ebm.obj import load_obj_age
    import pandas as pd
    from tqdm import tqdm
    import json

    folder = '../out/post_dyage'

    obj = load_obj_age(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0
    )

    with open(f'{folder}/Sim_Baseline.json', 'r') as f:
        y0s = json.load(f)

    mss = list()

    for i, state in tqdm(enumerate(y0s)):
        y0, pars = state['Y0'], state['Pars']
        y0 = np.array(y0)

        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)

        intvs = {
            'Baseline': compose_intv(p),
            'HighPPM': compose_intv(p, tx='HighPPM'),
        }

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'{folder}/Sim_IntvRel.csv')
