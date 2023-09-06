__author__ = 'Chu-Chang Ku'


exo = {
    'drt_act': 0
}

if __name__ == '__main__':
    from sim.ebm.obj import load_obj_baseline
    import pandas as pd
    import numpy as np
    from tqdm import tqdm

    suffix = 'free'

    obj = load_obj_baseline(
        folder_input=f'../data',
        file_prior='../data/prior.txt',
        year0=2000,
        exo=exo,
        suffix=suffix
    )

    post = pd.read_csv('../out/post_free/Post.csv')

    post = [dict(row) for i, row in post.iterrows()]

    mss = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)
        _, ms, _ = obj.Model.simulate_to_fit(p, t_eval=np.linspace(2018, 2040, 2040 - 2018 + 1))
        mss.append(ms.assign(Key=i))
    mss = pd.concat(mss)
    mss.to_csv('../out/post_free/Sim.csv')
