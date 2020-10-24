#
# dsdobjects/objectio.py
#   - copy and/or modify together with tests/test_objectio.py
#
import logging
log = logging.getLogger(__name__)

from .singleton import SingletonError
from .iupac_utils import reverse_wc_complement
from .dsdparser import (parse_seesaw_string, parse_seesaw_file,
                        parse_pil_string, parse_pil_file)

Domain = None
Complex = None
Reaction = None
Macrostate = None

def set_prototypes(): # Replace all objects with prototypes
    from .base_classes import DomainS as D
    from .base_classes import ComplexS as C
    from .base_classes import ReactionS as R
    from .base_classes import MacrostateS as M

    global Domain
    global Complex
    global Reaction
    global Macrostate

    Domain = D
    Complex = C
    Reaction = R
    Macrostate = M
    R.RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

class MissingObjectError(Exception):
    pass

class PilFormatError(Exception):
    pass

# ---- Load prototype objects ---- #
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

    domains = {}
    complexes = {}
    macrostates = {}
    det_reactions = []
    con_reactions = []
    for line in parsed_file :
        if ignore and line[0] in ignore:
            continue
        obj = read_pil_line(line)
        if Domain and isinstance(obj, Domain):
            domains[obj.name] = obj
            domains[(~obj).name] = ~obj
        elif Complex and isinstance(obj, Complex):
            complexes[obj.name] = obj
        elif Macrostate and isinstance(obj, Macrostate):
            macrostates[obj.name] = obj
        elif Reaction and isinstance(obj, Reaction) and obj.rtype == 'condensed':
            con_reactions.append(obj)
        elif Reaction and isinstance(obj, Reaction) and obj.rtype != 'condensed':
            det_reactions.append(obj)
    return domains, complexes, macrostates, det_reactions, con_reactions

def read_pil_line(raw):
    """ Interpret a single line of PIL input format.  """
    if isinstance(raw, str):
        [line] = parse_pil_string(raw)
    else:
        line = raw

    name = line[1]
    if line[0] == 'dl-domain':
        if Domain is None:
            raise MissingObjectError(f'Domain object not found: {Domain}')
        dlen = 5 if line[2] == 'short' else 15 if line[2] == 'long' else int(line[2]) 
        anon = Domain(name, length = dlen)
        return anon

    elif line[0] == 'sl-domain':
        if Domain is None:
            raise MissingObjectError(f'Domain object not found: {Domain}')
        if len(line) == 4:
            if int(line[3]) != len(line[2]):
                raise PilFormatError("Sequence/Length information inconsistent {line[3]} vs {len(line[2])}.")
        anon = Domain(name, length = len(line[2]))
        anon.sequence = line[2]
        comp = ~anon
        # Assume we are using DNA molecules and WC complements, for now.
        comp.sequence = reverse_wc_complement(line[2], material = 'DNA')
        return anon
 
    elif line[0] == 'kernel-complex':
        if Complex is None:
            raise MissingObjectError(f'Complex object not found: {Complex}')
        sequence, structure = resolve_kernel_loops(line[2])
        try : # to replace names with domain objects.
            sequence = [Domain(x) if x != '+' else '+' for x in sequence]
        except KeyError as err:
            raise PilFormatError(f"Cannot find domain: {err}.")
        cplx = Complex(sequence, structure, name = name)
        if len(line) > 3 :
            assert len(line[3]) == 3
            if hasattr(cplx, 'concentration') and cplx.concentration is not None:
                log.warning("Updating concentration for complex '{name}' to {line[3]}.")
            cplx.concentration = (line[3][0], float(line[3][1]), line[3][2])
        return cplx

    elif line[0] == 'resting-macrostate':
        if Macrostate is None:
            raise MissingObjectError('Macrostate object not found: {Macrostate}')
        try: # to replace names with complex objects.
            cplxs = [Complex(None, None, x) for x in line[2]]
        except KeyError as err:
            raise PilFormatError(f"Cannot find complex: {err}.")
        return Macrostate(complexes = cplxs, name = cplxs[0].name)

    elif line[0] == 'reaction':
        if Reaction is None:
            raise MissingObjectError(f'Reaction object not found: {Reaction}')
        reactants, products, rtype, rate, units, r = read_reaction(line)
        if rtype == 'condensed' :
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

        anon.set_rate_parameters(rate, units = units)
        #if anon.rateunits != units:
        #    raise PilFormatError(f"Rate units must be given in {anon.rateunits}, not: {units}.")
        return anon
    raise PilFormatError('Unknown keyword: {line[0]}')

