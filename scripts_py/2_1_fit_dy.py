__author__ = 'Chu-Chang Ku'

seed = 1167

year0 = 2000
n_iter = 300
n_round = 35
n_collect = 300

exo = {
    'drt_act': 0,
    'drt_trans': 0,
    'k_relapse_adj': 1,
    'rr_relapse_pub': 1,
    'r_acquire_dr': 0.03,
    # 'rr_beta_dr': 0.35,
}


if __name__ == '__main__':
    from sims_pars.fit import ApproxBayesComSMC
    import numpy.random as rd
    import os
    from sim.dy.obj import load_obj

    rd.seed(seed)

    obj = load_obj(
        folder_input='../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        exo=exo,
        suffix='cas_cdx',
        agp=True
    )

    # Fitting
    alg = ApproxBayesComSMC(n_iter=n_iter, max_round=n_round, parallel=True)
    alg.fit(obj)

    # Collect posterior
    post = alg.sample_posteriors(n_collect)

    # Output
    out_path = '../out/dyage'
    os.makedirs(out_path, exist_ok=True)

    post.Notes['Trace'].to_csv(f'{out_path}/Trace.csv')
    post.to_pred_df().to_csv(f'{out_path}/Pred.csv')
    post.to_df().to_csv(f'{out_path}/Post.csv')
