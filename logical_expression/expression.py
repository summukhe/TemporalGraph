import numpy as np
from copy import deepcopy
from .expression_utility import *


__version__ = "1.0"
__all__ = ['Numeral', 'Operator', 'Variable']


class Entity:
    def __init__(self, type, name):
        assert type in {'Operator', 'Numeral', 'Variable'}
        self.__type = type
        self.__name = name

    @property
    def type(self):
        return self.__type

    @property
    def name(self):
        return self.__name

    def eval(self):
        return 0

    @property
    def get_attached_variable(self):
        return []

    def __str__(self):
        return '%s' % self.name

    @property
    def depth(self):
        return 1

    @property
    def n_variables(self):
        return 1

    @property
    def n_operators(self):
        return 0


class Operator(Entity):
    def __init__(self, name, nvar, fn):
        assert type(fn).__name__ == 'function'
        Entity.__init__(self, type="Operator", name=name)
        assert nvar in {'unary', 'binary'}
        self.__type = nvar
        self.__fn = fn
        self.__children = [None] if type == 'unary' else [None, None]

    @property
    def depth(self):
        if not self.is_set:
            return 0

        if self.__type == 'unary':
            if isinstance(self.__children[0], Operator):
                return 1 + self.__children[0].depth
            else:
                return 1
        else:
            if isinstance(self.__children[0], Operator):
                ldepth = 1 + self.__children[0].depth
            else:
                ldepth = 1
            if isinstance(self.__children[1], Operator):
                rdepth = 1 + self.__children[1].depth
            else:
                rdepth = 1
            return np.maximum(ldepth, rdepth)

    @property
    def n_variables(self):
        return len(set([ent.name for ent in self.get_attached_variable]))

    @property
    def n_operators(self):
        if not self.is_set:
            return 1
        if self.__type == 'unary':
            return 1 + self.__children[0].n_operators
        else:
            n_l = self.__children[0].n_operators
            n_r = self.__children[1].n_operators
            return 1 + n_l + n_r

    @property
    def n_children(self):
        return 1 if self.__type == 'unary' else 2

    def __getitem__(self, item):
        return self.__children[item]

    def __setitem__(self, item, value):
        assert isinstance(value, Entity)
        self.__children[item] = value

    @property
    def is_set(self):
        if self.__type == 'unary':
            return self.__children[0] is not None
        else:
            return (self.__children[0] is not None) and (self.__children[1] is not None)

    @property
    def value(self):
        return self.eval()

    def eval(self):
        if not self.is_set:
            raise Exception('Operator not set')

        if self.__type == 'unary':
            return self.__fn(self.__children[0].value)
        else:
            return self.__fn(self.__children[0].value, self.__children[1].value)

    @property
    def get_attached_variable(self):
        data = []
        if self.is_set:
            if self.__type == 'unary':
                data = self.__children[0].get_attached_variable
            else:
                data = self.__children[0].get_attached_variable + self.__children[1].get_attached_variable
        return data

    def __str__(self):
        if self.__type == 'unary':
            return '(%s %s)' % (self.name, self.__children[0])
        else:
            return '(%s %s %s)' % (self.__children[0], self.name, self.__children[1])


class Numeral(Entity):
    def __init__(self, name, value):
        Entity.__init__(self, type="Numeral", name=name)
        assert is_supported_dtype(value)
        self.__value = deepcopy(value)

    @property
    def value(self):
        return self.eval()

    def eval(self):
        return deepcopy(self.value)

    @property
    def dtype(self):
        return get_dtype(self.__value)

    @property
    def shape(self):
        return shape(self.__value)

    def __str__(self):
        return '%s' % self.__value


class Variable(Entity):
    def __init__(self, name, variable):
        assert isinstance(name, str)
        Entity.__init__(self, type='Variable', name=name)
        self.__value = variable

    @property
    def value(self):
        return self.eval()

    @property
    def get_dtype(self):
        return get_dtype(self.__value)

    def eval(self):
        return self.__value[0]

    @property
    def get_attached_variable(self):
        return [self]

    def __str__(self):
        return '%s' % self.name

