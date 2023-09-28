from sim.healthcare import *


seed = 1167
year0 = 2000
poc_txi = 0.95

ss = {
    'Baseline': None,
    # 'TSwab': lambda pp: get_intv_tswab(pp),
    # 'POC': lambda pp: get_intv_poc(pp, target=0.8, p_txi_poc=poc_txi),
    # 'POC_Hi': lambda pp: get_intv_poc(pp, target=0.95, p_txi_poc=None)
}


def extract_cascade(sys, pars):
    p_ent = sys['p_ent']
    p_itt = pars['p_itt1']

    p_dx = sys['p_dx']

    p_txi = sys['p_txi']

    cas_itt = (p_ent * p_itt).sum()
    cas_dx = (p_ent * p_itt * p_dx).sum()
    cas_txi = (p_ent * p_itt * p_dx * p_txi).sum()

    p_d = p_itt * p_dx * p_txi

    gap_itt, gap_dx, gap_txi = 1 - cas_itt, 1 - cas_dx / cas_itt, 1 - cas_txi / cas_dx

    n_vis = 1 / cas_txi
    ks = 2.7 / n_vis
    p_d2 = p_d / ks

    return {
        'Key': i,
        'CasITT': cas_itt, 'CasDx': cas_dx, 'CasTxi': cas_txi,
        'GapITT': gap_itt, 'GapDx': gap_dx, 'GapTxi': gap_txi,
        'N_Vis': n_vis,
        'PDx_Pub': p_d[0], 'PDx_Eng': p_d[1], 'PDx_Pri': p_d[2],
        'PDx_Pub2': p_d2[0], 'PDx_Eng2': p_d2[1], 'PDx_Pri2': p_d2[2],
        'CasITT2': cas_itt / ks, 'CasDx2': cas_dx / ks, 'CasTxi2': cas_txi / ks,
    }


if __name__ == '__main__':
    import numpy as np
    from sim.inputs import load_inputs
    import pandas as pd

    inp = load_inputs('../pars')

    post = pd.read_csv('../out/post_dyage/Post.csv')

    cs = list()

    for i, pars in post.iterrows():
        pars = dict(pars)
        pars = inp.Cascade.prepare_pars(pars)

        for k, v in ss.items():
            d = extract_cascade(pars, pars)
            d['Scenario'] = k
            cs.append(d)

    cs = pd.DataFrame(cs)
    cs.to_csv('../out/post_dy/Sim_IntvDx_Cas2.csv')
