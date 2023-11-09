__author__ = 'Chu-Chang Ku'


year0 = 2000

if __name__ == '__main__':
    from sim.dy.obj import load_obj
    import pandas as pd
    import numpy as np
    import json
    from tqdm import tqdm

    obj = load_obj(
        folder_input='../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        suffix='cas_cdx',
        agp=True
    )

    folder = '../out/dyage'

    post = pd.read_csv(f'{folder}/Post.csv')
    post = [dict(row) for i, row in post.iterrows()]

    mss = list()
    particles = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)
        _, ms, _ = obj.Model.simulate_to_fit(p, t_eval=np.linspace(2015, 2040, 2040 - 2015 + 1))
        mss.append(ms.assign(Key=i))

        ys, _, _ = obj.Model.simulate_to_baseline(p)
        y0 = ys.y[:, -1]

        particles.append({
            'Y0': y0.tolist(),
            'Pars': pars,
            'Year0': year0,
            'Key': i
        })

    mss = pd.concat(mss)
    mss.to_csv(f'{folder}/Sim.csv')

    with open(f'{folder}/Sim_Baseline.json',  'w') as f:
        json.dump(particles, f)
