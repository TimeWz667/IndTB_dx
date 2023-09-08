from collections import namedtuple
from sim.inputs.cascade import *
from sim.inputs.demography import *
import pandas as pd
import pickle as pkl

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs']


Inputs = namedtuple("Inputs", ('Demography', 'Cascade', 'Targets'))


def load_inputs(root, agp=False, cs_suffix='free'):
    if not agp:
        with open(f'{root}/ind_all_70to35.pkl', 'rb') as f:
            src = pkl.load(f)
    else:
        with open(f'{root}/ind_{agp}_70to35.pkl', 'rb') as f:
            src = pkl.load(f)
    demo = Demography(src)

    cr = CasRepo.load(f'{root}/pars_cs_{cs_suffix}.json')
    tar = pd.read_csv(f'{root}/targets.csv')

    return Inputs(demo, cr, tar)
