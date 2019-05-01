import numpy as np
from .signal import *
from copy import deepcopy
from temporal_graph.network_analysis import WeightedGraph

__author__ = "Sumanta Mukherjee"
__all__ = ["DiffusionOnStaticGraph"]


class DiffusionOnStaticGraph:
    def __init__(self, dissipation_rate=0.):
        assert (dissipation_rate >= 0) and (dissipation_rate <= 1.)
        self.__eps = 1e-5
        self.__G = None
        self.__time = 0
        self.__signal = SignalInput()
        self.__vertex_history = {}
        self.__signal_nodes = []
        self.__active_nodes = set()
        self.__time_array = []
        self.__visit = dict()
        self.__dissipation_rate = dissipation_rate
        self.__max_degree = 1
        self.__blocked = set([])

    def start(self, g, start_nodes, pulse, blocked=[]):
        assert isinstance(g, WeightedGraph)
        assert isinstance(pulse, SignalInput)
        assert isinstance(blocked, list)
        self.__G = deepcopy(g)
        self.__time = 0
        self.__signal = pulse
        self.__vertex_history = {v: [0.] for v in self.__G.vertices()}
        self.__max_degree = np.max([self.__G.out_degree(v) for v in self.__G.vertices()])
        assert isinstance(start_nodes, list)
        for n in start_nodes:
            assert n in self.__vertex_history
        self.__time_array = list([self.__time])
        self.__signal_nodes = deepcopy(start_nodes)
        self.__active_nodes = set()
        self.__visit = dict()
        self.__blocked = set([v for v in blocked if self.__G.is_vertex(v)])

    def step(self, dt=1.):
        assert (self.__G is not None)
        assert (len(self.__signal_nodes) > 0)
        assert dt > 0
        for n in self.__vertex_history.keys():
            self.__vertex_history[n].append(self.__vertex_history[n][-1])

        for n in self.__signal_nodes:
            if n not in self.__blocked:
                s = self.__signal.do_signal(self.__time)
                if s > 0:
                    self.__active_nodes.add(n)
                    self.__vertex_history[n][-1] += s

        remove_nodes, add_nodes = set(), set()
        for n in sorted(self.__active_nodes):
            v_connected = [v for v in self.__G.out_neighbors(n) if v not in self.__blocked]
            e_wts = np.array([self.__G.weight(n, v) for v in v_connected], dtype=np.float)
            if np.sum(e_wts) > 0:
                d = self.__G.out_degree(n)
                f = (self.__max_degree - d) * 1.0 * self.__dissipation_rate / self.__max_degree
                self.__vertex_history[n][-1] = (1 - f) * self.__vertex_history[n][-1]
                node_potential = self.__vertex_history[n][-1]
                diff_potentials = np.clip([self.__vertex_history[v][-1] - node_potential for v in v_connected],
                                          a_min=self.__eps,
                                          a_max=None)
                diff_potentials = diff_potentials / np.sum(diff_potentials)
                p_wts = e_wts / np.sum(e_wts)
                assert f >= 0
                for i, v in enumerate(v_connected):
                    s = p_wts[i] * diff_potentials[i] * node_potential
                    self.__vertex_history[v][-1] += s
                    self.__vertex_history[n][-1] -= s
                    add_nodes.add(v)
            if self.__vertex_history[n][-1] < self.__eps:
                remove_nodes.add(n)
        for a in add_nodes:
            self.__active_nodes.add(a)

        for a in self.__active_nodes:
            if a not in self.__visit:
                self.__visit[a] = self.__time + dt

        for r in remove_nodes:
            self.__active_nodes.remove(r)

        self.__time = self.__time + dt
        self.__time_array.append(self.__time)

    def block(self, v):
        if self.__G.is_vertex(v):
            self.__blocked.add(v)

    def unblock(self, v):
        if v in self.__blocked:
            self.__blocked.remove(v)

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

    def get_smooth_history(self, v, time_window=5):
        assert time_window > 0 and time_window % 2 == 1
        if (v in self.__vertex_history) and (len(self.__time_array) > time_window):
            result = {'time': [], 'value': []}
            n = len(self.__time_array)
            t_med = int(time_window // 2)
            for i, t in enumerate(range(t_med, n - t_med)):
                result['time'].append(self.__time_array[t])
                result['value'].append(np.mean([self.__vertex_history[v][i+p] for p in range(time_window)]))
            return result
        return None

    def hitting_time(self, v):
        if v in self.__visit:
            return self.__visit[v]
        return None

    def information_carried(self, v):
        s = 0
        if v in self.__vertex_history:
            for i in range(len(self.__time_array)):
                if i > 0:
                    s += 0.5 * (self.__time_array[i] - self.__time_array[i-1]) * \
                         (self.__vertex_history[v][i] + self.__vertex_history[v][i-1])
        return s

    def time_to_maximum(self, v):
        if v in self.__vertex_history:
            return self.__time_array[np.argmax(self.__vertex_history[v])]
        return None

    def signal_strength(self, v, t=0):
        if (v in self.__vertex_history) and (t > 0) and (t <= self.__time):
            if t == 0:
                return self.__vertex_history[v][0]
            idx = 0
            for it in self.__time_array:
                if it >= t:
                    break
                idx = idx + 1

            vn, tn = self.__vertex_history[v][idx], self.__time_array[idx]
            vp, tp = self.__vertex_history[v][idx-1], self.__time_array[idx-1]
            l = (t - tp)/(tn - tp)
            return (1 - l) * vp + l * vn
        return None





