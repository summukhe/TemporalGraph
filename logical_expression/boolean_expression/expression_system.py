import numpy as np
from .expression import *
from ..expression import *
from .parse_boolean_expression import *


__version__ = "1.0"
__all__ = ['BooleanExpressionSystem']


class BooleanExpressionSystem:
    def __init__(self):
        self.__vars = dict()
        self.__expression = list()

    def is_registered_variable(self, varname):
        return varname in self.__vars

    @property
    def register_variables(self):
        return list(self.__vars.keys())

    def register_expression(self, expr_str):
        assert isinstance(expr_str, str)
        fexpr = build_parse_tree(expr_str)
        assert fexpr.is_valid()
        variables = [r.root_name for r in fexpr.leaves]
        for var in variables:
            if not self.is_registered_variable(var):
                self.__vars[var] = np.array([0])

        def recursive_expression_build(ex):
            if ex is not None:
                assert isinstance(ex, BooleanParseTree) and ex.is_valid()
                if ex.root_type == 'Variable':
                    vname = ex.root_name
                    return Variable(name=vname, variable=self.__vars[vname])
                elif ex.root_type == 'Operator':
                    oname = ex.root_name
                    if oname == 'and':
                        op = AndOperator(var1=recursive_expression_build(ex.left_child),
                                         var2=recursive_expression_build(ex.right_child))
                    elif oname == 'or':
                        op = OrOperator(var1=recursive_expression_build(ex.left_child),
                                        var2=recursive_expression_build(ex.right_child))
                    else:
                        op = NotOperator(variable=recursive_expression_build(ex.right_child))
                    return op

        expr = recursive_expression_build(fexpr)
        found = False
        for e in self.__expression:
            if check_expression_equality(expr, e):
                found = True
                break
        if not found:
            self.__expression.append(expr)

    def __len__(self):
        return len(self.__expression)

    def __getitem__(self, item):
        if item in self.__vars:
            return self.__vars[item][0]
        else:
            return None

    def __setitem__(self, key, value):
        if key in self.__vars:
            self.__vars[key][0] = value
        else:
            self.__vars[key] = np.array([value])

    def get_expression(self, index):
        if (index >= 0) and (index < len(self.__expression)):
            return self.__expression[index]
        else:
            return None

    def eval(self):
        return [expr.eval() for expr in self.__expression]


