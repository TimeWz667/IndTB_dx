__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    import numpy as np
    from sim.dy.intervention import compose_intv
    from sim.dy.obj import load_obj_baseline
    import pandas as pd
    from tqdm import tqdm
    import json

    folder = '../out/dy'

    obj = load_obj_baseline(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        suffix='cas_cdx_alt'
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
        }

        r0 = pars['r_relapse_te']
        p0 = r0 / (1 + r0)

        for pr in np.linspace(0, 1, 11):
            r1 = r0 * (1 - pr)
            p1 =r1 / (1 + r1)
            intvs[f'RelRed_{pr}'] = compose_intv(p, tx=f'RelRed_{p1}')

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'{folder}/Sim_IntvTx_RelRed.csv')
