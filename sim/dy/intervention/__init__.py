from collections import namedtuple
from sim.dy.intervention.dx import *
from sim.dy.intervention.tx import *
from sim.dy.intervention.vaccine import *
from sim.dy.intervention.acf import *


__all__ = ['compose_intv']

Interventions = namedtuple('Interventions', ('Dx', 'Tx', 'Txs', 'Vac', 'ACF', 'PPM'))


def compose_intv(pars, dx='None', tx='None', vac='None', acf='0_0_0', ppm=0):
    return Interventions(Dx=get_intv_dx(pars, dx),
                         Tx=get_intv_tx(tx),
                         Txs=TxBPaLM(),
                         Vac=get_intv_vac(pars, vac),
                         ACF=get_intv_acf(acf),
                         PPM=get_intv_ppm(ppm))
