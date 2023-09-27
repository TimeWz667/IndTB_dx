from collections import namedtuple
from sim.ebm.intervention.dx import *
from sim.ebm.intervention.tx import *
from sim.ebm.intervention.vaccine import *
from sim.ebm.intervention.acf import *


Interventions = namedtuple('Interventions', ('Dx', 'Tx', 'Vac', 'ACF', 'PPM'))


def compose_intv(pars, dx='None', tx='None', vac='None', acf='0_0', ppm=0):
    return Interventions(Dx=get_intv_dx(pars, dx),
                         Tx=get_intv_tx(pars, tx),
                         Vac=get_intv_vac(pars, vac),
                         ACF=get_intv_acf(pars, acf),
                         PPM=get_intv_ppm(pars, ppm))