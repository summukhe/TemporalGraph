import numpy as np

__author__ = "Sumanta Mukherjee"
__all__ = ['FFNormalizer']


class FFNormalizer:
    def __init__(self, cutoff=0, potential=True):
        self.__cutoff = cutoff
        self.__potential = potential

    def __call__(self, potential):
        if self.__potential:
            value = -1 * potential
        else:
            value = potential
        return np.clip(value, a_min=self.__cutoff, a_max=None)