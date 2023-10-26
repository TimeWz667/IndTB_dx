import json
import pickle as pkl
from sim.healthcare import *
import numpy as np
import numpy.random as rd


n_collect = 2000
n_vis = 2.7

if __name__ == '__main__':
    pars_tx = json.load(open('../pars/pars_tx.json', 'r'))
    prev = pars_tx['prev']
    pars_tx = pars_tx['pars']

    for cdx in ['cdx', 'cdx_alt']:
        pars_dx = json.load(open(f'../pars/pars_{cdx}.json', 'r'))

        res = {
            'N_Collect': n_collect,
            'N_Visit': n_vis,
            'Particles': [],
            'Prev': prev
        }

        for _ in range(n_collect):
            p_dx = rd.choice(pars_dx)
            p_tx = rd.choice(pars_tx)

            particle = dict()

            particle.update({k: v for k, v in p_dx.items() if k.startswith('sens_') or k.startswith('spec_')})
            particle['p_scanty'] = p_dx['p_scanty']
            particle['dur_pub'] = 0.5
            particle['dur_pri'] = p_tx['dur_pri']
            particle['p_pri_on_pub'] = p_tx['p_pri_on_pub']
            particle['system'] = sys = get_system(pars_dx[0], pars_tx[0])

            particle['p_dx'] = np.array([sys.Public.seek_care(1, 0).TruePos, sys.Engaged.seek_care(1, 0).TruePos,
                                         sys.Private.seek_care(1, 0).TruePos])
            particle['p_txi'] = np.array([0.9, 0.85, 0.85])
            particle['p_ent'] = sys.Entry
            particle['p_itt'] = 1 / (n_vis * particle['p_dx'] * particle['p_txi'] * particle['p_ent']).sum()
            particle['txi'] = np.array([p_tx['txi_pub'], p_tx['txi_eng'], p_tx['txi_pri']]).sum()
            particle['ppv'] = np.array([p_dx['ppv_pub'], p_dx['ppv_eng'], p_dx['ppv_pri']])

            res['Particles'].append(particle)

        with open(f'../pars/pars_cas_{cdx}.pkl', 'wb') as f:
            pkl.dump(res, f)
