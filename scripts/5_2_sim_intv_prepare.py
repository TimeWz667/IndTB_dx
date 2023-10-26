__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
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

    post = pd.read_csv(f'{folder}/Post.csv')
    post = [dict(row) for i, row in post.iterrows()]

    particles = list()

    for i, pars in tqdm(enumerate(post)):
        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)

        ys, _, _ = obj.Model.simulate_to_baseline(p)
        y0 = ys.y[:, -1]

        particles.append({
            'Y0': y0.tolist(),
            'Pars': pars,
            'Year0': year0,
            'Key': i
        })

    with open(f'{folder}/Sim_Baseline.json',  'w') as f:
        json.dump(particles, f)
