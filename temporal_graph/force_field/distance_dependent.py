import os
import logging
import numpy as np
import pandas as pd
from temporal_graph.pdb_processor import *

__version__ = "1.0"
__all__ = ['DistanceDependentPotential']


class DistanceDependentPotential:
    class _PairPotential:
        def __init__(self):
            self.__logger = logging.getLogger(name="temporal_network.DistanceDependentPotential")
            curr_dir = os.path.dirname(os.path.abspath(__file__))
            potential_file = os.path.join(curr_dir, 'data', 'energy_potential.csv')
            if not os.path.isfile(potential_file):
                self.__logger.error('Can not find energy_potential.csv for distance dependent forcefield computation!')
                raise Exception("File not found: %s" % potential_file)
            data = pd.read_csv(potential_file, header=0, sep=",")
            fields = ['amino1', 'amino2', 'f4', 'f3', 'f2', 'f1', 'f0']
            for f in fields:
                if f not in data.keys():
                    self.__logger.error('Can not find desired field (%s) in the energy_potential.csv' % f)
                    raise Exception('Missing field (%s) in energy_potential.csv' % f)
            self.__potential = dict()
            for i in data.index:
                self.__logger.debug(data.loc[i,:])
                amino1 = get_amino(data.loc[i, 'amino1']).name(one_letter_code=True)
                amino2 = get_amino(data.loc[i, 'amino2']).name(one_letter_code=True)
                parameters = np.array(data.loc[i, ['f4', 'f3', 'f2', 'f1', 'f0']])
                self.__logger.debug('Reading index %d (%s, %s) : [%s]' % (i, amino1, amino2, parameters))
                if amino1 not in self.__potential:
                    self.__potential[amino1] = dict()
                self.__potential[amino1][amino2] = np.poly1d(parameters)

        def __call__(self, amino1, amino2, distance):
            if distance < 3.5 or distance > 15:
                return 0
            aa1 = get_amino(amino1).name(one_letter_code=True)
            aa2 = get_amino(amino2).name(one_letter_code=True)
            self.__logger.debug('Requested potential for (%s <-> %s) pair' % (aa1, aa2))
            if aa1 not in self.__potential:
                self.__logger.error('Missing potential parameter for (%s)' % aa1)
                raise Exception('Missing potential parameters (%s)' % aa1)
            if aa2 not in self.__potential[aa1]:
                self.__logger.error('Missing pair potential parameter (%s %s)' % (aa1, aa2))
                raise Exception('Missing pair potential parameter (%s %s)' % (aa1, aa2))
            f = self.__potential[aa1][aa2](distance)
            return 0 if (f < 0) else f

    __instance = None

    @staticmethod
    def get_instance():
        if DistanceDependentPotential.__instance is None:
            DistanceDependentPotential.__instance = DistanceDependentPotential._PairPotential()
        return DistanceDependentPotential.__instance


