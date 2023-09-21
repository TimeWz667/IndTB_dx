from sim.ebm.intervention import compose_intv

__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    from sim.ebm.obj import load_obj_age
    import pandas as pd
    from tqdm import tqdm

    obj = load_obj_age(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0
    )

    post = pd.read_csv('../out/post_dyage/Post.csv')
    post = [dict(row) for i, row in post.iterrows()]

    mss = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)

        ys, _, _ = obj.Model.simulate_to_baseline(p)
        y0 = ys.y[:, -1]

        intvs = {
            'Baseline': compose_intv(p),
            'PAN-TB': compose_intv(p, tx='PAN-TB'),
            'LA-INJ': compose_intv(p, tx='LA-INJ')
        }

        for intv_key, intv in intvs.items():
            p1 = dict(p)
            _, ms, _ = obj.Model.simulate_onward(y0, p1, intv)
            mss.append(ms.assign(Key=i, Scenario=intv_key))
    mss = pd.concat(mss)
    print(mss)
    mss.to_csv(f'../out/post_dyage/Sim_IntvTx.csv')
