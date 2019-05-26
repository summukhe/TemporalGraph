import numpy as np


__version__ = "1.0"
__all__ = ['supported_numeric_type', 'get_dtype', 'is_array',
           'shape', 'is_bool', 'is_float', 'is_integer', 'is_supported_dtype']


def supported_numeric_type():
    return ['bool', 'float', 'integer']


def is_array(x):
    return isinstance(x, np.ndarray)


def is_bool(x):
    return isinstance(x, np.bool) or (is_array(x) and (x.dtype == np.bool))


def is_integer(x):
    return isinstance(x, np.int) or (is_array(x) and (x.dtype == np.int))


def is_float(x):
    return isinstance(x, np.float) or (is_array(x) and (x.dtype == np.float))


def is_supported_dtype(x):
    return is_bool(x) or is_integer(x) or is_float(x)


def shape(x):
    assert is_supported_dtype(x)
    if is_array(x):
        return x.shape
    return 1,


def get_dtype(x):
    if not is_supported_dtype(x):
        return None
    elif is_float(x):
        return 'float'
    elif is_bool(x):
        return 'bool'
    elif is_integer(x):
        return 'integer'
