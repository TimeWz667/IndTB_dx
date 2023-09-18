from sim.ebm.intv import Intervention


__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000
poc_txi = 0.95

ss = {
    'Baseline': None,
    'TSwab': lambda pp: get_intv_tswab(pp),
    'POC': lambda pp: get_intv_poc(pp, target=0.8, p_txi_poc=poc_txi),
    'POC_Hi': lambda pp: get_intv_poc(pp, target=0.95, p_txi_poc=None)
}

if __name__ == '__main__':
    from sim.healthcare import *
    from sim.ebm.obj import load_obj_baseline
    import pandas as pd
    from tqdm import tqdm

    obj = load_obj_baseline(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0
    )

    post = pd.read_csv('../out/post_dy/Post.csv')
    post = [dict(row) for i, row in post.iterrows()][:50]

    mss = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p, extra_sys=ss)

        ys, _, _ = obj.Model.simulate_to_baseline(p)
        y0 = ys.y[:, -1]

        intvs = dict()

        for k, v in ss.items():
            intvs[k] = {
                'Dx': {
                    'System': f'sys_{k}',
                    'PrTxi': p[f'sys_{k}']['p_txi'],
                    'Year0': 2025,
                    'PreFlight': 2
                }
            }

        intvs = {k: Intervention.parse_obj(v) for k, v in intvs.items()}

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'../out/post_dy/Sim_IntvDx.csv')
