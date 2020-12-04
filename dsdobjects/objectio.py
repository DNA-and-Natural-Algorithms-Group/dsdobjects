#
# dsdobjects/objectio.py
#   - copy and/or modify together with tests/test_objectio.py
#
import logging
log = logging.getLogger(__name__)

import gc
from .singleton import SingletonError
from .iupac_utils import reverse_wc_complement
from .complex_utils import strand_table_to_sequence
from .dsdparser import (parse_seesaw_string, parse_seesaw_file,
                        parse_pil_string, parse_pil_file)
from .base_classes import DomainS, StrandS, ComplexS, MacrostateS, ReactionS

Domain = None
Strand = None
Complex = None
Macrostate = None
Reaction = None

def clear_io_objects():
    global Domain
    global Strand
    global Complex
    global Macrostate
    global Reaction

    Domain = None
    Strand = None
    Complex = None
    Macrostate = None
    Reaction = None
    gc.collect()
    return

def set_io_objects(D = None, S = None, C = None, M = None, R = None):
    global Domain
    global Strand
    global Complex
    global Macrostate
    global Reaction
    Domain = DomainS if D is None else D
    Strand = StrandS if S is None else S
    Complex = ComplexS if C is None else C
    Macrostate = MacrostateS if M is None else M
    Reaction = ReactionS if R is None else R
    gc.collect()
    return

class PilFormatError(Exception):
    pass

def read_reaction(line):
    """ Interpret the parser output for a reaction line.
    """
    rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None
    rate = float(line[1][1][0]) if line[1] != [] and line[1][1] != [] else None
    error = float(line[1][1][1]) if line[1] != [] and line[1][1] != [] and len(line[1][1]) == 2 else None
    units = line[1][2][0] if line[1] != [] and line[1][2] != [] else None

    r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
    if rate is None:
        log.warning(f"Ignoring input reaction without a rate: {r}")
        return None, None, None, None, None, None
    elif rtype is None or rtype not in Reaction.RTYPES:
        log.warning(f"Ignoring input reaction '{rtype}': {r}")
        return None, None, None, None, None, None
    else :
        r = "[{} = {:12g} {}] {}".format(rtype, rate, units, r)
    return line[2], line[3], rtype, rate, units, r

def resolve_kernel_loops(loop):
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
            se, ss = resolve_kernel_loops(dom)
            sequen.extend(se)
            struct.extend(ss)
            sequen.append(old + '*' if old[-1] != '*' else old[:-1])
            struct.append(')')
    return sequen, struct

def read_pil(data, is_file = False, ignore = None):
    """ Read PIL file format.
    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool, optional): True if data is a path to a file, False otherwise
        ignore (list, optional): A list of identifiers that should be ignored.
    """
    parsed_file = parse_pil_file(data) if is_file else parse_pil_string(data)

    out = {'domains': dict(),
           'strands': dict(),
           'complexes': dict(),
           'macrostates': dict(),
           'det_reactions': set(),
           'con_reactions': set(),
           'other': []}

    for line in parsed_file :
        if ignore and line[0] in ignore:
            continue
        obj = read_pil_line(line)
        if isinstance(obj, Domain):
            out['domains'][obj.name] = obj
            comp = ~obj
            if obj.sequence is not None and comp.sequence is None:
                comp.sequence = reverse_wc_complement(obj.sequence, material = 'DNA')
            out['domains'][comp.name] = comp
            del comp # so important
        elif isinstance(obj, Strand):
            out['strands'][obj.name] = obj
        elif isinstance(obj, Complex):
            out['complexes'][obj.name] = obj
        elif isinstance(obj, Macrostate):
            out['macrostates'][obj.name] = obj
        elif isinstance(obj, Reaction) and obj.rtype == 'condensed':
            out['con_reactions'].add(obj)
        elif isinstance(obj, Reaction) and obj.rtype != 'condensed':
            out['det_reactions'].add(obj)
        else:
            assert isinstance(obj, list)
            out['other'].append(obj)
        del obj # so important
    return out

def read_pil_line(raw):
    """ Interpret a single line of PIL input format.  """
    if isinstance(raw, str):
        [line] = parse_pil_string(raw)
    else:
        line = raw

    name = line[1]
    if line[0] == 'dl-domain' and Domain is not None:
        dlen = 5 if line[2] == 'short' else 15 if line[2] == 'long' else int(line[2]) 
        anon = Domain(name, length = dlen)
        return anon

    elif line[0] == 'sl-domain' and Domain is not None:
        if len(line) == 4:
            if int(line[3]) != len(line[2]):
                raise PilFormatError("Sequence/Length information inconsistent {line[3]} vs {len(line[2])}.")
        anon = Domain(name, length = len(line[2]))
        anon.sequence = line[2]
        return anon

    elif line[0] == 'composite-domain' and Strand is not None:
        # This could be a strand definition or a composite domain.
        sequence = [Domain(d) for d in line[2]]
        anon = Strand(sequence, name)
        return anon
 
    elif line[0] == 'strand-complex' and Complex is not None:
        st = [list(Strand(None, name = s).sequence) for s in line[2]]
        sequence = strand_table_to_sequence(st)
        structure = line[3].replace(' ','')
        anon = Complex(sequence, list(structure), name = name)
        return anon

    elif line[0] == 'kernel-complex' and Complex is not None:
        sequence, structure = resolve_kernel_loops(line[2])
        try: # to replace names with domain objects.
            sequence = [Domain(x) if x != '+' else '+' for x in sequence]
        except SingletonError as err:
            # The PIL language allows shortcuts to specify composite domains
            # and then a kernel string using those composite domains. 
            # The below code tries all sorts of possibilties on what the
            # character in a sequence could actually refer to before it
            # finally throws an exception.
            def comp(name):
                return name[:-1] if name[-1] == '*' else name + '*'
            for e, d in enumerate(sequence):
                if isinstance(d, Domain) or d == '+':
                    # Happens with composite domains, see next statement e+1
                    continue
                try:
                    sequence[e] = Domain(d) # The default case.
                except SingletonError:
                    try:
                        # The character refers to a composite domain.
                        subseq = list(Strand(None, name = d).sequence)
                    except SingletonError:
                        try:
                            # The character refers to the complement of a composite domain.
                            complement = list(Strand(None, name = comp(d)).sequence)
                            subseq = [~d for d in reversed(complement)]
                        except SingletonError:
                            raise PilFormatError(f"Cannot find domain: {d}.")
                    for i, sd in enumerate(subseq):
                        assert isinstance(sd, Domain)
                        if i == 0:
                            sequence[e] = sd
                        else:
                            sequence.insert(e+i, sd)
                            structure.insert(e+i, structure[e])

        cplx = Complex(sequence, structure, name = name)
        if len(line) > 3 :
            assert len(line[3]) == 3
            if cplx.concentration is not None:
                log.warning(f"Updating concentration for complex '{name}' to {line[3]}.")
            cplx.concentration = (line[3][0], float(line[3][1]), line[3][2])
        return cplx

    elif line[0] == 'resting-macrostate' and Macrostate is not None:
        try: # to replace names with complex objects.
            cplxs = [Complex(None, None, x) for x in line[2]]
        except KeyError as err:
            raise PilFormatError(f"Cannot find complex: {err}.")
        return Macrostate(complexes = cplxs, name = name)

    elif line[0] == 'reaction' and Reaction is not None:
        reactants, products, rtype, rate, units, r = read_reaction(line)
        if rtype == 'condensed':
            try:
                reactants = [Macrostate(None, x) for x in reactants]
                products = [Macrostate(None, x) for x in products]
            except KeyError as err:
                raise PilFormatError(f"Cannot find resting complex: {err}.")
            anon = Reaction(reactants, products, rtype)
        else :
            try:
                reactants = [Complex(None, None, x) for x in reactants]
                products = [Complex(None, None, x) for x in products]
            except KeyError as err:
                raise PilFormatError(f"Cannot find complex: {err}.")
            anon = Reaction(reactants, products, rtype)

        anon.rate_constant = (rate, units)
        return anon
    else:
        log.warning(f'Cannot interpret line starting with {line[0]}')
        return line

