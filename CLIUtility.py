import os
import numpy as np

__version__ = "1.0"
__all__ = ['is_valid_file', 'is_integer']


def is_valid_file(parser, x):
    if not os.path.isfile(x):
        parser.error("File (%s) does not exists" % x)
    return x


def is_integer(parser, x):
    if isinstance(x, np.int) or int(x) != 0:
        return int(x)
    parser.error('Not an integer value (%s)' % x)


