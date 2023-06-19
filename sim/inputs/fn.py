import json
from collections import namedtuple
from sim.inputs.cascade import *
from sim.inputs.demography import *

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs']


Inputs = namedtuple("Inputs", ('Demography', 'Cascade'))


def load_inputs(root):
    with open(f'{root}/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    cr = CasRepo.load(f'{root}/pars_tx.json')

    return Inputs(demo, cr)
