__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    from sim.dy.obj import load_obj
    from sim.dy.intervention import compose_intv
    import pandas as pd
    import numpy as np
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

    for i, state in tqdm(enumerate(y0s[:5])):
        y0, pars = state['Y0'], state['Pars']
        y0 = np.array(y0)

        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)

        intvs = {
            'Baseline': compose_intv(p),
            # 'Dx_TSwab': compose_intv(p, dx='TSwab'),
            # 'Dx_POC': compose_intv(p, dx='POC'),
            # 'Dx_POC_Hi': compose_intv(p, dx='POC_Hi'),
            # 'Dx_TSwab_ITT': compose_intv(p, dx='TSwab_ITT'),
            # 'Dx_POC_ITT': compose_intv(p, dx='POC_ITT'),
            # 'Dx_POC_Hi_ITT': compose_intv(p, dx='POC_Hi_ITT'),
            # 'DxPPM_TSwab': compose_intv(p, dx='TSwab', ppm=0.9),
            # 'DxPPM_POC': compose_intv(p, dx='POC', ppm=0.9),
            # 'DxPPM_POC_Hi': compose_intv(p, dx='POC_Hi', ppm=0.9),
            'Tx_PAN-TB': compose_intv(p, tx='PAN-TB'),
            'Tx_LA-INJ': compose_intv(p, tx='LA-INJ'),
            # 'TxPPM_PAN-TB': compose_intv(p, tx='PAN-TB', ppm=0.9),
            # 'TxPPM_LA-INJ': compose_intv(p, tx='LA-INJ', ppm=0.9),
            # 'Vac_BCG': compose_intv(p, vac='BCG'),
            # 'Vac_M72': compose_intv(p, vac='M72'),
            # 'Vac_BCG-M72': compose_intv(p, vac='BCG-M72'),
            # 'Vac_Recurrence': compose_intv(p, vac='Recurrence'),
            # 'Mass_Xray_10_HRZE': compose_intv(p, acf='0.11_0.96_0'),
            # 'Mass_NAAT_20_HRZE': compose_intv(p, acf='0.22_0.86_0'),
            # 'Mass_Xray_10_PAN': compose_intv(p, acf='0.11_0.96_1', tx='Null_PAN-TB'),
            # 'Mass_NAAT_20_PAN': compose_intv(p, acf='0.22_0.86_1', tx='Null_PAN-TB'),
            # 'Combine_Lo': compose_intv(p, dx='TSwab', acf='0.11_0.96_0', vac='BCG'),
            # 'Combine_Hi': compose_intv(p, dx='TSwab', tx='LA-INJ', acf='0.22_0.86_0', vac='BCG-M72'),
            # 'CombinePPM_Lo': compose_intv(p, dx='TSwab', acf='0.11_0.96_0', vac='BCG', ppm=0.9),
            # 'CombinePPM_Hi': compose_intv(p, dx='TSwab', tx='LA-INJ', acf='0.22_0.86_0', vac='BCG-M72', ppm=0.9),
        }

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'{folder}/Sim_IntvAll0.csv')
