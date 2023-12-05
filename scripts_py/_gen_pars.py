__author__ = 'Chu-Chang Ku'


seed = 1167
year0 = 2000


if __name__ == '__main__':
    from sim.dy.obj import load_obj
    import numpy as np
    from tqdm import tqdm
    import json

    folder = '../out/dyage'

    obj = load_obj(
        folder_input='../pars',
        file_prior='../data/prior.txt',
        file_targets='../data/targets.csv',
        year0=year0,
        suffix='cas_cdx',
        agp=True
    )

    with open(f'{folder}/Sim_Baseline.json', 'r') as f:
        y0s = json.load(f)

    pss = list()

    for i, state in tqdm(enumerate(y0s[:100])):
        y0, pars = state['Y0'], state['Pars']
        y0 = np.array(y0)

        p = obj.serve(pars)
        p = obj.Cas.prepare_pars(p)
        pr = {k: v for k, v in p.items() if k in ['sens_s', 'sens_x', 'sens_x_sn', 'spec_s', 'spec_x', 'spec_x_sn',
                                                  'sens_cdx', 'spec_cdx', 'sens_cdx_bn', 'spec_cdx_bn',
                                                  'p_loss_sample',
                                                  'dur_pub', 'dur_pri', 'p_pri_on_pub',
                                                  'r_die_sym', 'rr_die_asym', 'r_die_tx',
                                                  'beta', 'rr_beta_dr', 'rr_sus_slat',
                                                  'p_primary', 'p_clear', 'r_clear', 'r_sc', 'r_react', 'r_relapse', 'r_cure',
                                                  'irr_25', 'irr_35', 'irr_45', 'irr_55', 'irr_65', 'rt_cs',
                                                  'r_acquire_dr', 'r_lat', 'p_cure', 'p_cure_dr',
                                                  't0_decline', 'r_act', 'r_relapse_te', 'r_onset', 'r_csi', 'r_recsi'
                                                  ]}

        pr['p_txi'] = p['p_txi'].tolist()

        pss.append({
            'y0': y0.tolist(),
            'pars': pr
        })

    with open('../R_ver/pars.json', 'w') as f:
        json.dump(pss, f)
