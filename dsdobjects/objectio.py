#
# dsdobjects/objectio.py
#   - copy and/or modify together with tests/test_objectio.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# Distributed under the MIT License, use at your own risk.
#
from __future__ import absolute_import, division, print_function

import logging
log = logging.getLogger(__name__)

# Memorymanagement
from dsdobjects.core import clear_memory
from dsdobjects.core import DSDObjectsError, DSDDuplicationError
from dsdobjects.core import DL_Domain, SL_Domain 
from dsdobjects.core import DSD_Complex, DSD_Reaction, DSD_Macrostate, DSD_StrandOrder

# Parsing and utils
from dsdobjects.utils import resolve_loops
from dsdobjects.parser import parse_seesaw_string, parse_seesaw_file
from dsdobjects.parser import parse_pil_string, parse_pil_file, PilFormatError

LogicDomain = None
Domain = None
Complex = None
Reaction = None
Macrostate = None

def set_prototypes(): # Replace all objects with prototypes
    from dsdobjects.prototypes import LogicDomain as LD
    from dsdobjects.prototypes import Domain as D
    from dsdobjects.prototypes import Complex as C
    from dsdobjects.prototypes import Reaction as R
    from dsdobjects.prototypes import Macrostate as M

    global LogicDomain
    global Domain
    global Complex
    global Reaction
    global Macrostate

    LogicDomain = LD
    Domain = D
    Complex = C
    Reaction = R
    Macrostate = M

class MissingObjectError(Exception):
    pass

# ---- Load prototype objects ---- #
def read_reaction(line):
    rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None
    rate = float(line[1][1][0]) if line[1] != [] and line[1][1] != [] else None
    error = float(line[1][1][1]) if line[1] != [] and line[1][1] != [] and len(line[1][1]) == 2 else None
    units = line[1][2][0] if line[1] != [] and line[1][2] != [] else None

    if rate is None:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        log.warning("Ignoring input reaction without a rate: {}".format(r))
        return None, None, None, None, None, None
    elif rtype is None or rtype not in Reaction.RTYPES:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        log.warning("Ignoring input reaction of with rtype='{}': {}".format(rtype, r))
        return None, None, None, None, None, None
    else :
        r = "[{} = {:12g} {}] {} -> {}".format(
                rtype, rate, units, ' + '.join(line[2]), ' + '.join(line[3]))

    return line[2], line[3], rtype, rate, units, r

def read_pil(data, is_file = False, ignore = None):
    """ Read PIL file format.

    Use dsdobjects parser to extract information. Load kinda.objects.

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
    """
    parsed_file = parse_pil_file(data) if is_file else parse_pil_string(data)

    dl_domains = {}
    sl_domains = {}
    complexes = {}
    macrostates = {}
    con_reactions = []
    det_reactions = []
    for line in parsed_file :
        if ignore and line[0] in ignore:
            continue
        obj = read_pil_line(line)
        if LogicDomain and isinstance(obj, LogicDomain):
            dl_domains[obj.name] = obj
            dl_domains[(~obj).name] = ~obj
        elif Domain and isinstance(obj, Domain):
            sl_domains[obj.name] = obj
            sl_domains[(~obj).name] = ~obj
        elif Complex and isinstance(obj, Complex):
            complexes[obj.name] = obj
        elif Macrostate and isinstance(obj, Macrostate):
            macrostates[obj.name] = obj
        elif Reaction and isinstance(obj, Reaction) and obj.rtype == 'condensed':
            con_reactions.append(obj)
        elif Reaction and isinstance(obj, Reaction) and obj.rtype != 'condensed':
            det_reactions.append(obj)
    domains = sl_domains if len(sl_domains) >= len(dl_domains) else dl_domains

    return domains, complexes, macrostates, det_reactions, con_reactions

def read_pil_line(raw):
    """ Read a single line of PIL input format.

    For example: Add a reaction to your system.

    """
    if isinstance(raw, str):
        line = parse_pil_string(raw)
        assert len(line) == 1
        line = line[0]
    else:
        line = raw

    name = line[1]
    if line[0] == 'dl-domain':
        if LogicDomain is None:
            raise MissingObjectError('No LogicDomain object found: {}'.format(LogicDomain))
        if line[2] == 'short':
            (dtype, dlen) = ('short', None)
        elif line[2] == 'long':
            (dtype, dlen) = ('long', None)
        else :
            (dtype, dlen) = (None, int(line[2]))

        anon = LogicDomain(name, dtype = dtype, length = dlen)
        comp = ~anon
        return anon

    elif line[0] == 'sl-domain':
        if Domain is None:
            raise MissingObjectError('No Domain object found: {}'.format(Domain))
        if len(line) == 4:
            if int(line[3]) != len(line[2]):
                raise PilFormatError(
                        "Sequence/Length information inconsistent {} vs ().".format(
                            line[3], len(line[2])))

        dtype = LogicDomain(name, length = len(line[2]))
        anon = Domain(dtype, line[2])
        comp = ~anon
        return anon
 
    elif line[0] == 'kernel-complex':
        if Complex is None:
            raise MissingObjectError('No Complex object found: {}'.format(Domain))
        sequence, structure = resolve_loops(line[2])
        DL_Domain.MEMORY['+'] = '+'
        SL_Domain.MEMORY['+'] = {'+': '+'}
        try : # to replace names with domain objects.
            sequence = list(map(lambda d : SL_Domain.MEMORY[d][d], sequence))
        except KeyError as err:
            try:
                sequence = list(map(lambda d : DL_Domain.MEMORY[d], sequence))
            except KeyError as err:
                raise PilFormatError("Cannot find domain: {}.".format(err))
        
        try:
            cplx = Complex(sequence, structure, name=name)
        except DSDDuplicationError as err:
            cplx = err.existing
            if cplx.name != name:
                raise DSDObjectsError("Complex {} exists under different name: {}.".format(name, cplx.name))

        if len(line) > 3 :
            assert len(line[3]) == 3
            if cplx.concentration is not None:
                log.warning("Updating concentration for {} to {}.".format(name, line[3]))
            cplx.concentration = (line[3][0], float(line[3][1]), line[3][2])
        return cplx


    elif line[0] == 'resting-macrostate':
        if Macrostate is None:
            raise MissingObjectError('No Macrostate object found: {}'.format(Macrostate))
        try: # to replace names with complex objects.
            cplxs = list(map(lambda c : 
                DSD_Complex.MEMORY[DSD_Complex.NAMES[c]], line[2]))
        except KeyError as err:
            raise PilFormatError("Cannot find complex: {}.".format(err))
        
        try:
            return Macrostate(name = name, complexes = cplxs)
        except DSDDuplicationError as err:
            return err.existing

    elif line[0] == 'reaction':
        if Reaction is None:
            raise MissingObjectError('No Reaction object found: {}'.format(Reaction))
        reactants, products, rtype, rate, units, r = read_reaction(line)

        if rtype == 'condensed' :
            try:
                reactants = list(map(lambda c : 
                    DSD_Macrostate.MEMORY[DSD_Macrostate.NAMES[c]], reactants))
                products  = list(map(lambda c : 
                    DSD_Macrostate.MEMORY[DSD_Macrostate.NAMES[c]], products))
            except KeyError as err:
                raise PilFormatError("Cannot find resting complex: {}.".format(err))
            anon = Reaction(reactants, products, rtype, rate)
        else :
            try:
                reactants = list(map(lambda c : 
                    DSD_Complex.MEMORY[DSD_Complex.NAMES[c]], reactants))
                products  = list(map(lambda c : 
                    DSD_Complex.MEMORY[DSD_Complex.NAMES[c]], products))
            except KeyError as err:
                raise PilFormatError("Cannot find complex: {}.".format(err))
            anon = Reaction(reactants, products, rtype, rate)

        if anon.rateunits != units:
            raise SystemExit("Rate units must be given in {}, not: {}.".format(
                        reaction.rateunits, units))
        return anon

    raise PilFormatError('unknown keyword: {}'.format(line[0]))

