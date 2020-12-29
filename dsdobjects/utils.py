#
# dsdobjects/utils.py
#   - copy and/or modify together with tests/test_utils.py
#
import re

# DEPRECATED  START
import warnings
from .complex_utils import split_complex_pt
from .complex_utils import strand_table_to_sequence
from .complex_utils import make_strand_table
from .complex_utils import make_loop_index as mli
from .complex_utils import make_pair_table as mpt
from .complex_utils import pair_table_to_dot_bracket as ptdb

class UtilityError(Exception):
    pass

def split_complex(stab, ptab):
    warnings.warn('dsdobjects: deprecated function call split_complex')
    def h():
        for st, pt in split_complex_pt(stab, ptab):
            nseq = strand_table_to_sequence(st)
            nsst = pair_table_to_dot_bracket(pt)
            yield (nseq, nsst)
    return list(h())

def make_lol_sequence(s):
    warnings.warn('dsdobjects: deprecated function make_lol_sequence: use make_strand_table')
    return make_strand_table(s)

def make_pair_table(*args, **kwargs):
    warnings.warn('dsdobjects: deprecated import')
    return mpt(*args, **kwargs)

def pair_table_to_dot_bracket(*args, **kwargs):
    warnings.warn('dsdobjects: deprecated import')
    return ptdb(*args, **kwargs)

def make_loop_index(s):
    warnings.warn('dsdobjects: deprecated import')
    return mli(s)
# DEPRECATED  END

def flint(n):
    """ int or float: of a number. """
    try:
        return int(float(n)) if float(n) == int(float(n)) else float(n)
    except OverflowError:
        return n

def convert_units(val, unit_in, unit_out):
    conc = {'M': 1, 'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12}
    time = {'days': 86400, 'hours': 3600, 'min': 60,
            's': 1, 'ms': 1e-3, 'us': 1e-6, 'ns': 1e-9}
    if unit_in in conc:
        return flint(val*conc[unit_in]/conc[unit_out])
    elif unit_in in time:
        return flint(val*time[unit_in]/time[unit_out])
    else:
        raise ValueError(f'Unknown unit for conversion: {unit_in}')

