import numpy as np
from copy import deepcopy
from .expression import *
from ..expression import *


__version__ = "1.0"
__all__ = ['build_random_boolean_expression', 'search_expr', 'expr_complexity']


def build_random_split_tree(vlist):
    assert isinstance(vlist, list) and len(vlist) > 0
    if len(vlist) < 3:
        return vlist
    else:
        rsplit_point = np.random.choice(len(vlist))
        while rsplit_point < 1:
            rsplit_point = np.random.choice(len(vlist))
        lsplit = vlist[:rsplit_point]
        rsplit = vlist[rsplit_point:]
        return [build_random_split_tree(lsplit), build_random_split_tree(rsplit)]


def recursive_build_random(vlist):
    assert isinstance(vlist, list) and len(vlist) == 2
    recusrive_left = isinstance(vlist[0], list)
    recursive_right = isinstance(vlist[1], list)
    if recusrive_left and len(vlist[0]) == 2:
        left_child = recursive_build_random(vlist[0])
    elif recusrive_left:
        left_child = NotOperator(vlist[0][0]) if np.random.random() < 0.5 else vlist[0][0]
    else:
        left_child = vlist[0]

    if recursive_right and len(vlist[1]) == 2:
        right_child = recursive_build_random(vlist[1])
    elif recursive_right:
        right_child = NotOperator(vlist[1][0]) if np.random.random() < 0.5 else vlist[1][0]
    else:
        right_child = vlist[1]
    return AndOperator(left_child, right_child) if np.random.random() < 0.5 else OrOperator(left_child, right_child)


def build_random_boolean_expression(variable_list,
                                    replicate_control=0.8):
    assert isinstance(variable_list, list)
    for var in variable_list:
        assert isinstance(var, Variable)
    n = len(variable_list)
    vlist = []
    for i in range(n):
        repl = int(np.random.random() > replicate_control) + 1
        for j in range(repl):
            vlist.append(variable_list[i])
    np.random.shuffle(vlist)
    vlist = build_random_split_tree(vlist)
    if len(vlist) == 1:
        return  vlist[0] if np.random.random() < 0.5 else NotOperator(vlist[0])
    else:
        assert len(vlist) == 2
        return recursive_build_random(vlist)


def expr_complexity(expr):
    return np.maximum(expr.n_variables, expr.n_operators)


def search_expr(y_values,
                x_values,
                dim_limit=3,
                nexpr=20,
                replicate_control=0.8,
                niter=100):
    assert isinstance(y_values, list)
    assert isinstance(x_values, list)
    assert len(x_values) == len(y_values)
    assert len(x_values) > 1
    assert dim_limit > 0
    assert nexpr > 1
    n_data = len(y_values)
    vnames = None
    for i in range(n_data):
        assert isinstance(x_values[i], dict)
        if vnames is None:
            vnames = set(x_values[i].keys())
        else:
            for k in x_values[i].keys():
                assert k in vnames

    dim_limit = np.minimum(dim_limit, len(vnames))
    variables, values = dict(), dict()
    for k in vnames:
        values[k] = np.array([0])
        variables[k] = Variable(name=k, variable=values[k])

    vnames = list(vnames)

    def score_expr(expr):
        score = 0
        for i in range(n_data):
            for v in vnames:
                values[v][0] = x_values[i][v]
            score = score + (y_values[i] != expr.eval())
        return score + expr_complexity(expr)

    def search_nexpr(nexpr):
        exprs, scores = list(), list()
        for i in range(nexpr):
            np.random.shuffle(vnames)
            nconsider = np.random.randint(low=1, high=dim_limit+1)
            vconsider = vnames[:nconsider]
            vlist = [variables[k] for k in vconsider]
            exprs.append(build_random_boolean_expression(variable_list=vlist, replicate_control=replicate_control))
        for expr in exprs:
            scores.append(score_expr(expr))
        score_args = np.argsort(scores)
        return [exprs[i] for i in score_args], [scores[i] for i in score_args]

    iter, stop = 0, False
    exprs , scores = search_nexpr(nexpr)

    while not stop:
        stop = iter > niter
        new_expr, new_scores = search_nexpr(nexpr)
        for i in range(nexpr):
            if scores[i] > new_scores[i]:
                scores[i] = new_scores[i]
                exprs[i] = new_expr[i]
                stop = False
        iter = iter + 1
    return exprs, scores


