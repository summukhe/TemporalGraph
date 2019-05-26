import numpy as np
from ..expression import *

__version__ = "1.0"
__all__ = ['AndOperator', 'OrOperator', 'NotOperator']


class AndOperator(Operator):
    def __init__(self, var1, var2):
        assert isinstance(var1, Variable) or isinstance(var1, Operator)
        assert isinstance(var2, Variable) or isinstance(var2, Operator)
        Operator.__init__(self, name='and', nvar='binary', fn=lambda x, y: np.minimum(x, y))
        self.__setitem__(0, var1)
        self.__setitem__(1, var2)


class OrOperator(Operator):
    def __init__(self, var1, var2):
        assert isinstance(var1, Variable) or isinstance(var1, Operator)
        assert isinstance(var2, Variable) or isinstance(var2, Operator)
        Operator.__init__(self, name='or', nvar='binary', fn=lambda x, y: np.maximum(x, y))
        self.__setitem__(0, var1)
        self.__setitem__(1, var2)


class NotOperator(Operator):
    def __init__(self, variable):
        assert isinstance(variable, Variable) or isinstance(variable, Operator)
        Operator.__init__(self, name='not', nvar='unary', fn=lambda x: -1 * x)
        self.__setitem__(0, variable)



