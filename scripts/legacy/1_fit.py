__author__ = 'Chu-Chang Ku'


n_iter = 300
n_round = 30
n_collect = 300
seed = 1167


exo = {
    'drt_act': 0,
    'drt_trans': 0,
    'rt_cs': 0.015
}


if __name__ == '__main__':
    from sims_pars.fit import ApproxBayesComSMC
    from sim.ebm.obj import load_obj_baseline
    import os

    suffix = 'free'

    obj = load_obj_baseline(
        folder_input=f'../data',
        file_prior='../data/prior.txt',
        year0=2000,
        exo=exo,
        suffix=suffix
    )

    # Fitting
    alg = ApproxBayesComSMC(n_iter=n_iter, max_round=n_round, parallel=True)
    alg.fit(obj)

    # Collect posterior
    post = alg.sample_posteriors(n_collect)
    print(post)

    # Output
    out_path = f'../out/post_{suffix}'
    os.makedirs(out_path, exist_ok=True)

    post.Notes['Trace'].to_csv(f'{out_path}/Trace.csv')
    post.to_pred_df().to_csv(f'{out_path}/Pred.csv')
    post.to_df().to_csv(f'{out_path}/Post.csv')
