__author__ = 'Chu-Chang Ku'


year0 = 2000

if __name__ == '__main__':
    from sim.ebm.obj import load_obj_baseline
    import pandas as pd
    from tqdm import tqdm
    from dy.intervention import compose_intv

    obj = load_obj_baseline(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        suffix='cas_cdx'
    )

    out_folder = '../out/post_dy_c0'

    post = pd.read_csv(f'{out_folder}/Post.csv')

    post = [dict(row) for i, row in post.iterrows()][:150]

    mss1 = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)
        ys, _, _ = obj.Model.simulate_to_baseline(p)

        y0 = ys.y[:, -1]

        intvs = {
            'Baseline': compose_intv(p),
            'Dx_TSwab': compose_intv(p, dx='TSwab'),
            'Dx_POC': compose_intv(p, dx='POC'),
            'Dx_POC_Hi': compose_intv(p, dx='POC_Hi'),
        }
        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss1.append(ms.assign(Key=i, Scenario=intv_key))

    mss1 = pd.concat(mss1)
    mss1.to_csv(f'{out_folder}/Sim_Intv.csv')
