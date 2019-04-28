from temporal_graph.network_analysis import WeightedGraph
from .signal import *
from copy import deepcopy
import numpy as np

__author__ = "Sumanta Mukherjee"
__all__ = ["InformationDiffusionOnStaticGraph"]


class InformationDiffusionOnStaticGraph():
    def __init__(self):
        self.__G = None
        self.__time = 0
        self.__signal = SignalInput()
        self.__vertex_history = {}
        self.__signal_nodes = []
        self.__active_nodes = set()
        self.__time_array = []

    def start(self, g, start_nodes, pulse):
        assert isinstance(g, WeightedGraph)
        assert isinstance(pulse, SignalInput)
        self.__G = deepcopy(g)
        self.__time = 0
        self.__signal = pulse
        self.__vertex_history = {v: [0.] for v in self.__G.vertices()}
        assert isinstance(start_nodes, list)
        for n in start_nodes:
            assert n in self.__vertex_history
        self.__time_array = list([self.__time])
        self.__signal_nodes = deepcopy(start_nodes)
        self.__active_nodes = set()

    def step(self, dt=1.):
        assert (self.__G is not None)
        assert (len(self.__signal_nodes) > 0)
        assert dt > 0
        for n in self.__vertex_history.keys():
            self.__vertex_history[n].append(self.__vertex_history[n][-1])

        for n in self.__signal_nodes:
            s = self.__signal.do_signal(self.__time)
            if s > 0:
                self.__active_nodes.add(n)
                self.__vertex_history[n][-1] += s

        remove_nodes, add_nodes = set(), set()
        for n in self.__active_nodes:
            v_connected = self.__G.out_neighbors(n)
            e_wts = np.array([self.__G.weight(n, v) for v in v_connected], dtype=np.float)
            e_wts = e_wts / np.sum(e_wts)
            for i, v in enumerate(v_connected):
                self.__vertex_history[v][-1] += self.__vertex_history[n][-1] * e_wts[i]
                add_nodes.add(v)
            remove_nodes.add(n)
        for a in add_nodes:
            self.__active_nodes.add(a)
        for r in remove_nodes:
            self.__active_nodes.remove(r)

        self.__time = self.__time + dt
        self.__time_array.append(self.__time)

    def curr_time(self):
        return self.__time

    def vertices(self):
        if self.__G is not None:
            return self.__G.vertices()
        return []

    def edges(self):
        if self.__G is not None:
            return self.__G.edges()
        return []

    def get_history(self, v):
        if v in self.__vertex_history:
            return {'time': self.__time_array,
                    'value': self.__vertex_history[v]}
        return None


