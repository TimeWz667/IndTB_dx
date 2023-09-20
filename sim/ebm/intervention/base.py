from abc import ABCMeta, abstractmethod
from collections import namedtuple
import numpy as np
from scipy.interpolate.interpolate import interp1d

__author__ = 'Chu-Chang Ku'
__all__ = ['IntvSingleton', 'Interventions']


Interventions = namedtuple('Interventions', ('Dx', 'Tx', 'Vac', 'ACF'))


class IntvSingleton(metaclass=ABCMeta):
    def __init__(self, pars):
        pass

    @abstractmethod
    def modify(self, t, y, pars):
        pass


if __name__ == '__main__':
    import numnpy as np

    x = np.linspace(0,)


    fn = interp1d(
        x = np.linspace()
    )
