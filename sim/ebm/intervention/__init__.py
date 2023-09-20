from collections import namedtuple
from sim.ebm.intervention.dx import *
from sim.ebm.intervention.tx import *
# from sim.ebm.intervention.vaccine import *


Interventions = namedtuple('Interventions', ('Dx', 'Tx', 'Vac', 'ACF'))


def compose_intv(pars, dx='None', tx='None', vac='None'):
    return Interventions(Dx=get_intv_dx(pars, dx),
                         Tx=get_intv_tx(pars, tx), Vac=None, ACF=None)
