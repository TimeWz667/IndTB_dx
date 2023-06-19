from abc import ABCMeta, abstractmethod

__author__ = 'Chu-Chang Ku'
__all__ = ['Process']


class Process(metaclass=ABCMeta):
    def __init__(self, pid, keys):
        self.PID = pid
        self.Keys = keys

    def find_dya(self, t, ya, pars, calc, **kwargs):
        y, aux = ya
        self.calculate_calc(t, y, pars, calc, **kwargs)
        return self.compose_dya(ya, calc), calc

    @abstractmethod
    def compose_dya(self, ya, calc: dict):
        pass

    def measure(self, mea: dict, t, ya, pars, calc: dict, **kwargs):
        pass

    @abstractmethod
    def calculate_calc(self, t, y, pars, calc, **kwargs):
        pass
