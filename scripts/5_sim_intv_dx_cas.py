from sim.healthcare import *


seed = 1167
year0 = 2000
poc_txi = 0.95

ss = {
    'Baseline': None,
    'TSwab': lambda pp: get_intv_tswab(pp),
    'POC': lambda pp: get_intv_poc(pp, target=0.8, p_txi_poc=poc_txi),
    'POC_Hi': lambda pp: get_intv_poc(pp, target=0.95, p_txi_poc=None)
}


def extract_cascade(sys, pars):
    p_ent = sys['sys'].Entry
    p_itt = pars['p_itt1']

    p_dx = sys['sys'].seek_care(1, 0)
    p_dx = np.array([d.TruePos for d in p_dx.values()]) / p_ent

    p_txi = sys['p_txi']

    cas_itt = (p_ent * p_itt).sum()
    cas_dx = (p_ent * p_itt * p_dx).sum()
    cas_txi = (p_ent * p_itt * p_dx * p_txi).sum()

    gap_itt, gap_dx, gap_txi = 1 - cas_itt, 1 - cas_dx / cas_itt, 1 - cas_txi / cas_dx

    return {
        'Key': i,
        'CasITT': cas_itt, 'CasDx': cas_dx, 'CasTxi': cas_txi,
        'GapITT': gap_itt, 'GapDx': gap_dx, 'GapTxi': gap_txi
    }


if __name__ == '__main__':
    import numpy as np
    from sim.inputs import load_inputs
    import pandas as pd

    inp = load_inputs('../pars')

    post = pd.read_csv('../out/post_dy/Post.csv')

    cs = list()

    for i, pars in post.iterrows():
        pars = dict(pars)
        pars = inp.Cascade.prepare_pars(pars, extra_sys=ss)

        for k, v in ss.items():
            d = extract_cascade(pars[f'sys_{k}'], pars)
            d['Scenario'] = k
            cs.append(d)

    cs = pd.DataFrame(cs)
    cs.to_csv('../out/post_dy/Sim_IntvDx_Cas.csv')
