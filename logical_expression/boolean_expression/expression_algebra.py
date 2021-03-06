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
                nexpr=10,
                replicate_control=0.8,
                niter=100,
                return_entropy=False):
    assert isinstance(y_values, list)
    assert isinstance(x_values, list)
    assert len(x_values) == len(y_values)
    assert len(x_values) > 1
    assert dim_limit > 0
    assert nexpr > 1
    n_data, vnames = len(y_values), None
    n_memory = np.maximum(100, nexpr)

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

    def sort_score_ordered(exprs, scores):
        assert isinstance(exprs, list)
        assert isinstance(scores, list)
        assert len(exprs) == len(scores)
        score_args = np.argsort(scores)
        return [exprs[i] for i in score_args], [scores[i] for i in score_args]

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
        return sort_score_ordered(exprs, scores)

    iter, stop = 0, False
    exprs , scores = search_nexpr(n_memory)

    while not stop:
        stop = iter > niter
        new_expr, new_scores = search_nexpr(n_memory)
        for i in range(nexpr):
            if scores[i] > new_scores[i]:
                found = False
                for expr in exprs:
                    if check_expression_equality(new_expr[i], expr):
                        found = True
                if not found:
                    scores[i] = new_scores[i]
                    exprs[i] = new_expr[i]
                    stop = False
        iter = iter + 1
    mod_scores = [s - expr_complexity(exprs[i]) for i, s in enumerate(scores)]
    exprs, scores = sort_score_ordered(exprs, mod_scores)
    exprs = exprs[:nexpr]
    scores = scores[:nexpr]
    coverage = [0 for i in range(n_data)]
    entropy = [[] for x in exprs]
    for i in range(n_data):
        for v in vnames:
            values[v][0] = x_values[i][v]
        for j, expr in enumerate(exprs):
            if expr.eval() == y_values[i]:
                coverage[i] += 1
                entropy[j].append(y_values[i])

    if return_entropy is False:
        return exprs, scores, coverage

    def calc_entropy(x):
        assert isinstance(x, list)
        n = len(x)
        count = dict()
        for v in x:
            if v not in count:
                count[v] = 0
            count[v] += 1./n
        e = 0
        for v in count:
            e = e - count[v]*np.log(count[v])
        return e

    return exprs, scores, coverage, [calc_entropy(x) for x in entropy]

