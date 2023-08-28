import json
from collections import namedtuple
from sim.inputs.cascade import *
from sim.inputs.demography import *

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs']


Inputs = namedtuple("Inputs", ('Demography', 'Cascade', 'Targets'))


def load_inputs(root):
    with open(f'{root}/pars_demo.json', 'r') as f:
        demo = Demography(json.load(f))
    cr = CasRepo.load(f'{root}/pars_cs.json')
    tar = cr.Targets

    return Inputs(demo, cr, tar)
