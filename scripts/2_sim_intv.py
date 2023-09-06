from sim.ebm.intv import Intervention

__author__ = 'Chu-Chang Ku'


intvs = {
    'Baseline': {},
    'RedRel': {'DeRel': {'Target': 0.02}},
    'MassScreen': {'ACF': {'Coverage': 0.1053605, 'Sens': 0.86}},
    'Vac': {
        'Vac': {
            'Efficacy': 0.5,
            'Coverage': 0.5,
            'Year0': 2025
        }
    },
    'TSwab': {'Swab': {'Uptake': 1, 'XpertAccess': 0.9}},
    'Txi': {'TxIni': {'Target': 0.98}},
    'All': {
        'DeRel': {'Target': 0.02},
        'ACF': {'Coverage': 0.1053605, 'Sens': 0.86},
        'Vac': {
                    'Efficacy': 0.5,
                    'Coverage': 0.5,
                    'Year0': 2025
                },
        'Swab': {'Uptake': 1, 'XpertAccess': 0.9},
        'TxIni': {'Target': 0.98}
    }
}
intvs = {k: Intervention.parse_obj(v) for k, v in intvs.items()}

if __name__ == '__main__':
    from sim.ebm.obj import load_obj_baseline
    import pandas as pd
    from tqdm import tqdm

    suffix = 'free'

    obj = load_obj_baseline(
        folder_input=f'../data',
        file_prior='../data/prior.txt',
        year0=2000,
        exo={},
        suffix=suffix
    )

    post = pd.read_csv('../out/post_free/Post.csv')

    post = [dict(row) for i, row in post.iterrows()]

    mss = list()

    for i, pars in tqdm(enumerate(post)):
        pars['beta'] *= 1.13
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)
        ys, _, _ = obj.Model.simulate_to_baseline(p)
        y0 = ys.y[:, -1]

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'../out/post_{suffix}/Sim_Intv.csv')
