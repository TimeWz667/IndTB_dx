from abc import ABCMeta, abstractmethod

__author__ = 'Chu-Chang Ku'
__all__ = ['Process']


class Process(metaclass=ABCMeta):
    def __init__(self, pid, keys):
        self.PID = pid
        self.Keys = keys

    @abstractmethod
    def find_dya(self, t, y, da, pars, calc, intv=None):
        pass
