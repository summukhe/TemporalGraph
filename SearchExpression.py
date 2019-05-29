import os
import sys
sys.path.append( os.path.dirname( os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from logical_expression import *


if __name__ == "__main__":
    filename = 'data/data_extract_updated.csv'
    x = pd.read_csv(filename, header=0, sep=',')
    all_hdr = list(x.keys())
    assert 'rna' in all_hdr
    var_cols = [hdr for hdr in all_hdr if hdr not in {'resid', 'mutation', 'atp', 'rna', 'helicase'}]
    x_values, y_values = list(), list()

    def process_data(v):
        if v > 1:
            return 1
        elif v < 0:
            return -1
        else:
            return 0

    for i in x.index:
        y_values.append(x.loc[i, 'rna'])
        x_values.append({v: process_data(x.loc[i, v]) for v in var_cols})

    exprs, scores, coverage = search_expr(y_values=y_values,
                                          x_values=x_values,
                                          dim_limit=3,
                                          nexpr=20,
                                          niter=500,
                                          replicate_control=1.0)

    print(coverage)
    for i, expr in enumerate(exprs):
        print('%s (%d) => %.1f' % (expr, expr_complexity(expr), 100 * (1 - (scores[i]/len(y_values)))))

