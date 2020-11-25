#
# dsdobjects/complex_utils.py
#   - copy and/or modify together with tests/test_complex_utils.py
#
import warnings
from functools import reduce
from itertools import chain, groupby

class SecondaryStructureError(Exception):
    pass

def make_pair_table(ss, strand_break = '+', ignore = set('.')):
    """ Return a secondary structure in form of pair table:

    Args:
      ss (str or list): secondary structure in dot-bracket format.
      cut (str, optional): character that defines a cut-point. Defaults to '+'.
      ignore (set, optional): a list of characters that are ignored. Defaults to ['.']

    Example:
                      0, 0  0, 1  0, 2   0, 3    0, 4      1, 0    1, 1
      "...((+))" -> [[None, None, None, (1, 1), (1, 0)], [(0, 4), (0, 3)]]

    Raises:
       SecondaryStructureError: Too many closing parenthesis ')' in secondary structure.
       SecondaryStructureError: Too many opening parenthesis '(' in secondary structure.
       SecondaryStructureError: Unexpected character in sequence: "{}".

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
        if char == strand_break:
            strand_index += 1
            domain_index = 0
            strand = []
            pair_table.append(strand)
            continue
        if char == '(':
            strand.append(None)
            stack.append((strand_index, domain_index))
            domain_index += 1
        elif char == ')':
            try:
                loc = stack.pop()
            except IndexError as e:
                raise SecondaryStructureError("Too few closing parenthesis ')' in secondary structure.")
            strand.append(loc)
            pair_table[loc[0]][loc[1]] = (strand_index, domain_index)
            domain_index += 1
        elif char in set(ignore): # unpaired
            strand.append(None)
            domain_index += 1
        else:
            raise SecondaryStructureError(f"Unexpected character in sequence: '{char}'.")
    if len(stack) > 0:
        raise SecondaryStructureError("Too few opening parenthesis '(' in secondary structure.")
    return pair_table 

def pair_table_to_dot_bracket(pt, strand_break = '+', join = False):
    """ Inverse of the make_pair_table function.
    """
    assert len(strand_break) == 1
    out = ''
    for si, strand in enumerate(pt):
        if out: 
            out += strand_break
        for di, pair in enumerate(strand):
            if pair is None:
                out += '.'
            else:
                locus = (si, di)
                if locus < pair:
                    out += '('
                else:
                    out += ')'
    return out if join else list(out)

def make_strand_table(seq, strand_break = '+'):
    """ Returns a sequence in a format corresponding to pair tables or loop indices. """
    assert len(strand_break) == 1
    if isinstance(seq, list):
        return [list(g) for k, g in groupby(seq, key = lambda x: x != strand_break) if k]
    else:
        return [list(s) for s in seq.split(strand_break)]

def strand_table_to_sequence(st, strand_break = '+', join = False):
    """ Returns a sequence in a format corresponding to pair tables or loop indices. 

    Note, the join argument is actually important. It is natural to represent 
    nucleotide sequences as string, but for domain sequences that is not very
    useful. 
    """
    if join:
        return f'{strand_break}'.join(''.join(s) for s in st)
    else:
        return reduce(lambda a, b: a + [strand_break] + b, st)

def make_loop_index(pt, components = False):
    """ Return the loop index of a secondary structure.

    Each nucleotide gets a number assigned. Two nucleotides with the
    same number are either paired or able to pair to each other, i.e.
    they are in the same loop.

    Note, the parameter components changes the output format. If it is
    True, then disconnected complexes are allowed and the second output
    gives for each strand a pair of [exterior loop start, exterior loop stop].
    This can be used to split complexes into their connected components later.

    Returns:
        * A list of lists with the loop index for every nucleotide.
        * A set of loop indices that correspond to exterior loops.
    """
    loop_index = []
    exterior = set()
    myext = []

    stack = []
    (cl, nl) = (0, 0)
    for si, strand in enumerate(pt):
        loop = []
        loop_index.append(loop)
        ext = [cl, None]
        for di, pair in enumerate(strand):
            loc = (si, di)
            if pair is None:
                pass
            elif loc < pair: # '('
                nl += 1
                cl = nl
                stack.append(loc)
            loop.append(cl)
            if pair and pair < loc: # ')'
                _ = stack.pop()
                try:
                    ploc = stack[-1]
                    cl = loop_index[ploc[0]][ploc[1]]
                except IndexError:
                    cl = 0
        ext[1] = cl
        myext.append(ext)
        if cl in exterior:
            if components is False:
                raise SecondaryStructureError('Complexes not connected.')
        else:
            exterior.add(cl) 
    return (loop_index, exterior) if not components else (loop_index, myext)

def split_complex_db(seq, sst, join = False):
    """Yields connected complexes from dot bracket format. """
    stab = make_strand_table(seq)
    ptab = make_pair_table(sst)
    for st, pt in split_complex_pt(stab, ptab):
        nseq = strand_table_to_sequence(st, join = join)
        nsst = pair_table_to_dot_bracket(pt, join = join)
        yield (nseq, nsst)
    return

def split_complex_pt(stab, ptab):
    """ Split a complex into unconnected components.
    
    Args:
        stab: A sequence in "strand table" format.
        ptab: A structure in "pair table" format.

    Yields:
        Unconnected components.
    """
    def splice(i,j):
        innerss = stab[i:j+1]
        innerpt = [list(map(lambda x: (x if x is None else (x[0] - i, 
                                                            x[1])), st)) for st in ptab[i:j+1]]
        outerss = stab[0:i] + stab[j+1:]
        outerpt = [list(map(lambda x: (x if x is None else (x[0] if x[0] < i else x[0] - (j+1-i), 
                                                            x[1])), st)) for st in ptab[0:i] + ptab[j+1:]]
        return (innerss, innerpt), (outerss, outerpt)

    li, ext = make_loop_index(ptab, components = True)
    seen = {0: 0}
    for j, (fr, to) in enumerate(ext):
        assert fr in seen
        assert seen[fr] == j
        if j == len(ext)-1:
            assert to in seen
            yield (stab, ptab)
            break
        if to in seen:
            i = seen[to]
            (iss, ipt), (oss, opt) = splice(i, j)
            for (s, p) in chain(split_complex_pt(iss, ipt), 
                                split_complex_pt(oss, opt)):
                yield (s, p)
            break
        seen[to] = j + 1
    return

def wrap(x, m):
    """ Wraps x so that 0 <= wrap(x, m) < m. (x can be negative.) """
    return (x % m + m) % m

def rotate_complex_db(seq, sst, turns = None, join = False):
    # Would this be much faster when using strings?
    stab = make_strand_table(seq)
    ptab = make_pair_table(sst)
    assert all(len(x) == len(y) for x, y in zip(stab, ptab))
    for (st, pt) in rotate_complex_pt(stab, ptab, turns = turns):
        nseq = strand_table_to_sequence(st, join = join)
        nsst = pair_table_to_dot_bracket(pt, join = join)
        yield (nseq, nsst)
    return

def rotate_complex_pt(stab, ptab, turns = None):
    """ Returns all complex rotations starting with the current one.
    
    Note, you can set turns = 1 for a single forced rotation.
    """
    if turns is None:
        turns = len(ptab)

    def rotate_locus(x, n):
        return x if x is None else (wrap(x[0]+n, len(ptab)), x[1])

    if turns > 0:
        if len(ptab) > 1 and turns != len(ptab):
            stab = [stab[-1]] + stab[:-1]
            ptab = [[rotate_locus(x, 1) for x in y] for y in [ptab[-1]] + ptab[:-1]]
        yield (stab, ptab)
        for (s, p) in rotate_complex_pt(stab, ptab, turns - 1):
            yield (s, p)
    return

def rotate_complex_once(seq, sst, turns = None):
    # This function turns out to be much faster than the other two...
    stack = []
    if "+" in seq:
        p = seq.index('+')
        seq = seq[p + 1:] + ["+"] + seq[:p]
        nstr = list(sst)
        stack = []
        for i in range(p):
            if nstr[i] == "(":
                stack.append(i)
            elif nstr[i] == ")":
                try:
                    stack.pop()
                except IndexError:
                    raise SecondaryStructureErrro("Unbalanced parenthesis.")
        for i in stack:
            nstr[i] = ")"
        stack = []
        for i in reversed(range(p + 1, len(nstr))):
            if nstr[i] == ")":
                stack.append(i)
            elif nstr[i] == "(":
                try :
                    stack.pop()
                except IndexError:
                    raise SecondaryStructureError("Unbalanced parenthesis.")
        for i in stack:
            nstr[i] = "("
        sst = nstr[p + 1:] + ["+"] + nstr[:p]
    return seq, sst

