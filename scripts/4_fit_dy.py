

n_iter = 50
n_round = 15
n_collect = 50
seed = 1167

year0 = 2000

exo = {
    'drt_act': 0,
    'drt_trans': 0
}


if __name__ == '__main__':
    from sims_pars.fit import ApproxBayesComSMC
    import numpy.random as rd
    import os
    from sim.ebm.obj import load_obj_baseline

    rd.seed(seed)

    obj = load_obj_baseline(
        folder_input=f'../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        exo=exo
    )

    # Fitting
    alg = ApproxBayesComSMC(n_iter=n_iter, max_round=n_round, parallel=True)
    alg.fit(obj)

    # Collect posterior
    post = alg.sample_posteriors(n_collect)

    # Output
    out_path = '../out/post_dy'
    os.makedirs(out_path, exist_ok=True)

    post.Notes['Trace'].to_csv(f'{out_path}/Trace.csv')
    post.to_pred_df().to_csv(f'{out_path}/Pred.csv')
    post.to_df().to_csv(f'{out_path}/Post.csv')
