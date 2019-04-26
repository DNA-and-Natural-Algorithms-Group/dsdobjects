#
# dsdobjects/utils.py
#   - copy and/or modify together with tests/test_utils.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# Distributed under the MIT License, use at your own risk.
#

from functools import reduce

class DSDUtilityError(Exception):
    """
    dnaobjects.basic error class.
    """

    def __init__(self, message, *kargs):
        if kargs:
            self.message = "{} [{}]".format(message, ', '.join(map(str,kargs)))
        else :
            self.message = message
        super(DSDUtilityError, self).__init__(self.message)

def make_pair_table(ss, strand_break='+', ignore=set('.')):
    """Return a secondary struture in form of pair table:

    Args:
      ss (str or list): secondary structure in dot-bracket format.
      cut (str, optional): character that defines a cut-point. Defaults to '+'.
      ignore (set, optional): a list of characters that are ignored. Defaults to ['.']

    Example:
                             0,0  0,1  0,2   0,3   0,4     1,0   1,1
             "...((+))" -> [[None,None,None,(1,0),(1,1)],[(0,4),(0,3)]]

    Raises:
       DSDUtilityError: Too many closing parenthesis ')' in secondary structure.
       DSDUtilityError: Too many opening parenthesis '(' in secondary structure.
       DSDUtilityError: Unexpected character in sequence: "{}".

    Returns:
      [list]: A pair-table as list of lists.
    """

    assert len(strand_break) == 1
    assert '.' in ignore

    # Return value
    pair_table = []

    # Helpers
    stack = []
    strand_index = 0
    domain_index = 0

    strand = []
    pair_table.append(strand)
    for char in ss:
        if char == strand_break :
            strand_index += 1
            domain_index = 0
            strand = []
            pair_table.append(strand)
            continue

        if char == '(' :
            strand.append(None)
            stack.append((strand_index, domain_index))
            domain_index += 1
 
        elif char == ")" :
            try:
                loc = stack.pop()
            except IndexError as e:
                raise DSDUtilityError("Too few closing parenthesis ')' in secondary structure.")
            strand.append(loc)
            pair_table[loc[0]][loc[1]] = (strand_index, domain_index)
            domain_index += 1

        elif char in set(ignore): # unpaired
            strand.append(None)
            domain_index += 1

        else :
            raise DSDUtilityError("Unexpected character in sequence: '{}'.".format(char))

    if len(stack) > 0 :
        raise DSDUtilityError("Too few opening parenthesis '(' in secondary structure.")

    return pair_table 

def pair_table_to_dot_bracket(pt, strand_break='+'):
    """
    Inverse of the make_pair_table function.
    """
    assert len(strand_break) == 1
    out = []
    for si, strand in enumerate(pt):
        if out: out.append(strand_break)
        for di, pair in enumerate(strand):
            if pair is None:
                out.append('.')
            else :
                locus = (si, di)
                if locus < pair :
                    out.append('(')
                else :
                    out.append(')')
    return out

def make_lol_sequence(seq):
    indices = [-1] + [i for i, x in enumerate(seq) if x == "+"]
    indices.append(len(seq))
    return [seq[indices[i - 1] + 1 : indices[i]] for i in range(1, len(indices))]

def make_loop_index(ptable):
    """
    number loops and assign each position its loop-number
    handy for checking which pairs can be added

    Returns:
        * A list of lists with the loop index for every nucleotide
        * A set of loop indices that correspond to exterior loops.
    """
    loop_index = []
    exterior = set()

    stack = []
    (cl, nl) = (0, 0)

    for si, strand in enumerate(ptable):
        loop = []
        loop_index.append(loop)
        for di, pair in enumerate(strand):
            loc = (si, di)
            if pair is None:
                pass
            elif loc < pair : # '('
                nl += 1
                cl = nl
                stack.append(loc)
            loop.append(cl)
            if pair and pair < loc : # ')'
                _ = stack.pop()
                try :
                    ploc = stack[-1]
                    cl = loop_index[ploc[0]][ploc[1]]
                except IndexError:
                    cl = 0
        # strand break
        if cl in exterior :
            raise DSDUtilityError('Complexes not connected.')
        else :
            exterior.add(cl)
    # ptable end
    return loop_index, exterior

def split_complex(lol_seq, ptable):
    """
    NOTE: modifies its arguments!
    """
    loop_index = []
    exterior = set()

    complexes = dict()
    splice = []

    stack = []
    (cl, nl) = (0, 0)

    # Identify exterior loops using the loop index function,
    # rewrite ptable to contain characters instead of paired locus.
    for si, strand in enumerate(ptable):
        loop = []
        loop_index.append(loop)
        for di, pair in enumerate(strand):
            loc = (si, di)
            if pair is None:
                ptable[si][di] = '.'
            elif loc < pair : # '('
                nl += 1
                cl = nl
                stack.append(loc)
                ptable[si][di] = '('
            loop.append(cl)
            if pair and pair < loc : # ')'
                _ = stack.pop()
                try :
                    ploc = stack[-1]
                    cl = loop_index[ploc[0]][ploc[1]]
                except IndexError:
                    cl = 0
                ptable[si][di] = ')'

        # store the strand_index where a given loop-cut starts
        if cl in complexes :
            start = complexes[cl]
            end = si
            splice.append((start,end))
            complexes[cl] = si
        else :
            complexes[cl] = si
    # for ptable end

    parts = []
    for (s,e) in splice:
        new_seq = lol_seq[s+1:e+1]
        new_pt = ptable[s+1:e+1]

        lol_seq[s+1:e+1] = [[] for x in range(s+1,e+1)]
        ptable[s+1:e+1] = [[] for x in range(s+1,e+1)]

        # change format from lol to regular list format
        new_seq=reduce(lambda a,b:a+['+']+b, new_seq)
        new_pt=reduce(lambda a,b:a+['+']+b, new_pt)
        parts.append((new_seq, new_pt))

    # append the remainder (or original) sequence
    new_seq = lol_seq[:]
    new_pt = ptable[:]
    new_seq=reduce(lambda a,b:a+['+']+b if b else a+b,new_seq)
    new_pt=reduce(lambda a,b:a+['+']+b if b else a+b,new_pt)
    parts.append((new_seq, new_pt))
    return parts

def resolve_loops(loop):
    """ Return a sequence, structure pair from kernel format.
    """
    sequen = []
    struct = []
    for dom in loop :
        if isinstance(dom, str):
            sequen.append(dom)
            if dom == '+' :
                struct.append('+')
            else :
                struct.append('.')
        elif isinstance(dom, list):
            struct[-1] = '('
            old = sequen[-1]
            se, ss = resolve_loops(dom)
            sequen.extend(se)
            struct.extend(ss)
            sequen.append(old + '*' if old[-1] != '*' else old[:-1])
            struct.append(')')
    return sequen, struct

def convert_units(val, unit_in, unit_out):
    conc = {'M':1, 'mM':1e-3, 'uM':1e-6, 'nM':1e-9, 'pM':1e-12}
    time = {'ns':1e-9, 'us':1e-6, 'ms':1e-3, 's':1, 'min':60, 'hours':3600, 'days':86400}
    if unit_in in conc:
        return val*conc[unit_in]/conc[unit_out]
    elif unit_in in time:
        return val*time[unit_in]/time[unit_out]
    else:
        raise DSDUtilityError('Unknown unit for conversion: {}'.format(unit_in))


