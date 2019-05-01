#
# dsdobjects/prototypes.py
#   - copy and/or modify together with tests/test_prototypes.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# Distributed under the MIT License, use at your own risk.
#
# Commonly useful DSD object definitions.
#   - functionality from here may be incorporated into base_classes if generally useful.
#   - you should be able to copy that file as is into your project if you need custom changes.
#   - please consider providing thoughts about missing functionality
#
from __future__ import absolute_import, division, print_function

from dsdobjects.core import clear_memory
from dsdobjects.core import DSDObjectsError, DSDDuplicationError
from dsdobjects.core import SequenceConstraint # just pass it on ...
from dsdobjects.core import DL_Domain, SL_Domain 
from dsdobjects.core import DSD_Complex, DSD_Reaction, DSD_Macrostate, DSD_StrandOrder

from dsdobjects.utils import split_complex, resolve_loops
from dsdobjects.parser import parse_pil_string, parse_pil_file, PilFormatError

class LogicDomain(DL_Domain):
    """
    Represents a single domain. We allow several options for specifying domain
    properties. Domains might have an explicit integer (bp) length, or may be
    designated as short or long. If the latter method is used, the code will use
    the relevant constant as the integer domain length.
    """

    def __new__(cls, name, dtype=None, length=None):
        # The new method returns the present instance of an object, if it exists
        self = DL_Domain.__new__(cls)
        try:
            super(LogicDomain, self).__init__(name, dtype, length)
        except DSDDuplicationError as e :
            other = e.existing
            if dtype and (other.dtype != dtype) :
                raise DSDObjectsError('Conflicting dtype assignments for {}: "{}" vs. "{}"'.format(
                    name, dtype, other.dtype))
            elif length and (other.length != length) :
                raise DSDObjectsError('Conflicting length assignments for {}: "{}" vs. "{}"'.format(
                    name, length, other.length))
            return e.existing
        self._nucleotides = None
        return self

    def __init__(self, name, dtype=None, length=None):
        # Remove default initialziation to get __new__ to work
        pass

    @property
    def nucleotides(self):
        return self._nucleotides

    @nucleotides.setter
    def nucleotides(self, value):
        self._nucleotides = value

    @property
    def complement(self):
        # If we initialize the complement, we need to know the class.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = LogicDomain(cname, self.dtype, self.length)
        return self._complement

    def can_pair(self, other):
        """
        Returns True if this domain is complementary to the argument.
        """
        return self == ~other

    @property
    def identity(self):
        """
        Returns the identity of this domain, which is its name without a
        complement specifier (i.e. A and A* both have identity A).
        """
        return self._name[:-1] if self._name[-1] == '*' else self._name

    @property
    def is_complement(self):
        """
        Returns true if this domain is a complement (e.g. A* rather than A),
        false otherwise.
        """
        return self._name[-1:] == '*'

class Domain(SL_Domain):
    def __init__(self, dtype, sequence, variant=''):
        super(Domain, self).__init__(dtype, sequence, variant='')
    
    @property
    def complement(self):
        dtcomp = self._dtype.complement
        if dtcomp.name not in SL_Domain.MEMORY:
            d = Domain(dtcomp, sequence = 'N' * len(dtcomp))
        if len(list(SL_Domain.MEMORY[dtcomp.name].values())) > 1:
            raise NotImplementedError('complementarity not properly implemented')
        return list(SL_Domain.MEMORY[dtcomp.name].values())[0]

class Complex(DSD_Complex):
    """
    Complex prototype object. 

    Overwrites some functions with new names, adds some convenient stuff..
    """

    PREFIX = 'e'

    @staticmethod
    def clear_memory(memory=True, names=True, ids=True):
        if memory:
            DSD_Complex.MEMORY = dict()
        if names:
            DSD_Complex.NAMES = dict()
        if ids:
            DSD_Complex.ID = dict()

    def __init__(self, sequence, structure, name='', prefix='', memorycheck=True):
        try :
            if not prefix :
                prefix = Complex.PREFIX
            super(Complex, self).__init__(sequence, structure, name, prefix, memorycheck)
        except DSDObjectsError :
            backup = 'enum' if prefix != 'enum' else 'proto'
            super(Complex, self).__init__(sequence, structure, name, backup, memorycheck)
            logging.warning('Complex name existed, prefix has been changed to: {}'.format(backup))
        
        self.concentration = None # e.g. (initial, 5, nM)
        assert self.is_domainlevel_complement

    def get_strand(self, loc):
        """
        Returns the strand at the given index in this complex
        """
        if loc is not None :
            return self._strands[loc]
        return None

    @property
    def available_domains(self):
        ad = []
        for (x,y) in self.exterior_domains:
            ad.append((self.get_domain((x,y)), x, y))
        return ad

    @property
    def pk_domains(self):
        pd = []
        for (x,y) in self.enclosed_domains:
            pd.append((self.get_domain((x,y)), x, y))
        return pd

    def split(self):
        """ Split Complex into disconneted components.
        """
        if self.is_connected:
            return [self]
        else :
            ps = self.lol_sequence
            pt = self.pair_table
            parts = split_complex(ps, pt)
            cplxs = []
            # assign new_complexes
            for (se,ss) in parts:
                try:
                    cplxs.append(Complex(se, ss))
                except DSDDuplicationError as e:
                    cplxs.append(e.existing)
            return sorted(cplxs)

class Macrostate(DSD_Macrostate):
    pass

class Reaction(DSD_Reaction):
    RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    def __init__(self, *kargs, **kwargs):
    #def __init__(self, reactants, products, rtype=None, rate=None, memorycheck=True):
        super(Reaction, self).__init__(*kargs, **kwargs)
        if self._rtype not in Reaction.RTYPES:
            try:
                del DSD_Reaction.MEMORY[self.canonical_form]
            except KeyError:
                pass
            raise DSDObjectsError('Reaction type not supported! ' + 
            'Set supported reaction types using Reaction.RTYPES')

    def full_string(self, molarity='M', time='s'):
        """Prints the reaction in PIL format.
        Reaction objects *always* specify rate in /M and /s.  """

        if self.rate :
            newunits = [molarity] * (self.arity[0] - 1) + [time]
            newrate = self.rateformat(newunits)
            rate = newrate.constant
            assert newunits == newrate.units
            units = ''.join(map('/{}'.format, newrate.units))
        else :
            rate = float('nan')
            units = ''

        if self.rtype :
            return '[{:14s} = {:12g} {:4s} ] {} -> {}'.format(self.rtype, rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))
        else :
            return '[{:12g} {:4s} ] {} -> {}'.format(rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

    def ptreact(self):
        """ 
        Find a common pairtable representation for input and output.

        Needs thorough testing!

        Note: It seems like we can only do that if either len(reactants) == 1 or
        len(products)==1. Only then we have sufficient constraints on the strand
        order. For example, a reaction of strands: VXY + Z -> ... -> YZ + VX we
        might chose VXYZ for the strand order, even though the intermediate has
        VXZY.

        Returns:
            StrandOrder, pairtable-reactants, pairtable-products.

        """
        so = None  # The common strand order.
        pt1 = None # Pair table of reactants
        pt2 = None # Pair table of products

        rotations = 0
        if len(self.reactants) == 1:
            cplx = self.reactants[0]
            if isinstance(cplx, Macrostate):
                cplx = cplx.canonical
            try:
                so = StrandOrder(cplx.sequence, prefix='so_')
            except DSDDuplicationError as e : 
                so = e.existing
                rotations = e.rotations
            
            if rotations:
                for e, rot in enumerate(cplx.rotate()):
                    if e == rotations:
                        pt1 = rot.pair_table
            else:
                pt1 = cplx.pair_table

        elif len(self.products) == 1:
            cplx = self.products[0]
            if isinstance(cplx, Macrostate):
                cplx = cplx.canonical
            try:
                so = StrandOrder(cplx.sequence, prefix='so_')
            except DSDDuplicationError as e : 
                so = e.existing
                rotations = e.rotations
 
            if rotations:
                for e, rot in enumerate(cplx.rotate()):
                    if e == rotations:
                        pt2 = rot.pair_table
            else:
                pt2 = cplx.pair_table

        else :
            raise NotImplementedError

        # So now that we got a valid StrandOrder, we need to represent the 
        # other side as a disconnected Complex within that StrandOrder.
        # complexes = get_complexes_from_other_side()
        cxs = self.reactants if pt2 else self.products

        if any(map(lambda c: isinstance(c, Macrostate), cxs)):
            cxs = list(map(lambda x: x.canonical, cxs))

        assert len(cxs) == 2 # or better <= ?

        for rot1 in cxs[0].rotate():
            for rot2 in cxs[1].rotate():
                rotations = None
                try:
                    so2 = StrandOrder(rot1.sequence + ['+'] + rot2.sequence)
                except DSDDuplicationError as e : 
                    so2 = e.existing
                    rotations = e.rotations
                if so == so2:
                    fake = Complex(rot1.sequence + ['+'] + rot2.sequence,
                                   rot1.structure + ['+'] + rot2.structure, 
                                   memorycheck=False)
                    if rotations:
                        for x in range(len(so)-rotations):
                            fake = fake.rotate_once()
                    if pt1: pt2 = fake.pair_table
                    else :  pt1 = fake.pair_table

        # What have we got?
        return (so, pt1, pt2)

class StrandOrder(DSD_StrandOrder):
    pass

# ---- Load prototype objects ---- #
def read_reaction(line):
    rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None
    rate = float(line[1][1][0]) if line[1] != [] and line[1][1] != [] else None
    error = float(line[1][1][1]) if line[1] != [] and line[1][1] != [] and len(line[1][1]) == 2 else None
    units = line[1][2][0] if line[1] != [] and line[1][2] != [] else None

    if rate is None:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        logging.warning("Ignoring input reaction without a rate: {}".format(r))
        return None, None, None, None, None, None
    elif rtype is None or rtype not in Reaction.RTYPES:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        logging.warning("Ignoring input reaction of with rtype='{}': {}".format(rtype, r))
        return None, None, None, None, None, None
    else :
        r = "[{} = {:12g} {}] {} -> {}".format(
                rtype, rate, units, ' + '.join(line[2]), ' + '.join(line[3]))

    return line[2], line[3], rtype, rate, units, r

def read_pil(data, is_file = False):
    """ Read PIL file format.

    Use dsdobjects parser to extract information. Load kinda.objects.

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
    """
    parsed_file = parse_pil_file(data) if is_file else parse_pil_string(data)

    dl_domains = {'+' : '+'} # saves some code
    sl_domains = {'+' : '+'} # saves some code

    def assgn_dl_domain(name, dt, dl):
        """ Initialize both the domain and its complement. """
        if name not in dl_domains:
            dl_domains[name] = LogicDomain(name, dtype = dt, length = dl)
            cname = name[:-1] if dl_domains[name].is_complement else name + '*'
            dl_domains[cname] = ~dl_domains[name]
        return dl_domains[name]

    def assgn_sl_domain(dtype, sequence, variant=None):
        """ Initialize both the domain and its complement. 
            Does not support variants at this point.
        """
        if variant is None:
            name = dtype.name
            assert name in dl_domains
        else: 
            raise NotImplementedError

        if name not in sl_domains:
            sl_domains[name] = Domain(dtype, sequence)
            cname = name[:-1] if dl_domains[name].is_complement else name + '*'
            sl_domains[cname] = ~sl_domains[name]
        return sl_domains[name]

    sequences = {} # strand order?
    complexes = {}
    macrostates = {}
    con_reactions = []
    det_reactions = []
    for line in parsed_file :
        name = line[1]
        if line[0] == 'dl-domain':
            if line[2] == 'short':
                (dtype, dlen) = ('short', None)
            elif line[2] == 'long':
                (dtype, dlen) = ('long', None)
            else :
                (dtype, dlen) = (None, int(line[2]))

            if name not in dl_domains:
                _ = assgn_dl_domain(name, dtype, dlen)

        elif line[0] == 'sl-domain':
            if len(line) == 4:
                if int(line[3]) != len(line[2]):
                    raise PilFormatError(
                            "Sequence/Length information inconsistent {} vs ().".format(
                                line[3], len(line[2])))

            if name not in sl_domains:
                dldom = assgn_dl_domain(name, dt=None, dl=len(line[2]))
                sldom = assgn_sl_domain(dldom, sequence = line[2])
            else:
                assert sl_domains[name].sequence == line[2]

        elif line[0] == 'composite-domain':
            pass

        elif line[0] == 'strand-complex':
            pass
 
        elif line[0] == 'kernel-complex':
            sequence, structure = resolve_loops(line[2])

            try : # to replace names with domain objects.
                sequence = list(map(lambda d : sl_domains[d], sequence))
            except KeyError as err:
                try:
                    sequence = list(map(lambda d : dl_domains[d], sequence))
                except KeyError as err:
                    raise PilFormatError("Cannot find domain: {}.".format(err))
            
            complexes[name] = Complex(sequence, structure, name=name)

            if len(line) > 3 :
                assert len(line[3]) == 3
                complexes[name]._concentration = tuple(line[3])

        elif line[0] == 'resting-macrostate':
            try: # to replace names with complex objects.
                cplxs = list(map(lambda c : complexes[c], line[2]))
            except KeyError as err:
                raise PilFormatError("Cannot find complex: {}.".format(err))
            macrostates[name] = Macrostate(name = name, complexes = cplxs)

        elif line[0] == 'reaction':
            reactants, products, rtype, rate, units, r = read_reaction(line)

            if rtype == 'condensed' :
                try:
                    reactants = list(map(lambda c : macrostates[c], line[2]))
                    products  = list(map(lambda c : macrostates[c], line[3]))
                except KeyError as err:
                    raise PilFormatError("Cannot find resting complex: {}.".format(err))
                rxn = Reaction(reactants, products, rtype, rate)
                con_reactions.append(rxn)
            else :
                try:
                    reactants = list(map(lambda c : complexes[c], reactants))
                    products  = list(map(lambda c : complexes[c], products))
                except KeyError as err:
                    raise PilFormatError("Cannot find complex: {}.".format(err))

                rxn = Reaction(reactants, products, rtype, rate)
                det_reactions.append(rxn)

            if rxn.rateunits != units:
                raise SystemExit("Rate units must be given in {}, not: {}.".format(reaction.rateunits, units))

        else :
            print('# Ignoring keyword: {}'.format(line[0]))

    domains = sl_domains if len(sl_domains) >= len(dl_domains) else dl_domains

    return domains, complexes, macrostates, det_reactions, con_reactions

def read_pil_line(raw):
    line = parse_pil_string(raw)

    assert len(line) == 1
    line = line[0]

    name = line[1]
    if line[0] == 'dl-domain':
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
        
        if len(line) > 3 :
            assert len(line[3]) == 3
            print("WARNING: Ignoring concentration.")

        try:
            return Complex(sequence, structure, name=name)
        except DSDDuplicationError as err:
            return err.existing


    elif line[0] == 'resting-macrostate':
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
        reactants, products, rtype, rate, units, r = read_reaction(line)

        if rtype == 'condensed' :
            try:
                reactants = list(map(lambda c : 
                    DSD_Macrostate.MEMORY[DSD_Macrostate.NAMES[c]], line[2]))
                products  = list(map(lambda c : 
                    DSD_Macrostate.MEMORY[DSD_Macrostate.NAMES[c]], line[3]))
            except KeyError as err:
                raise PilFormatError("Cannot find resting complex: {}.".format(err))
            anon = Reaction(reactants, products, rtype, rate)
        else :
            try:
                reactants = list(map(lambda c : complexes[c], 
                    DSD_Complex.MEMORY[DSD_Complex.NAMES[c]], reactants))
                products  = list(map(lambda c : complexes[c], 
                    DSD_Complex.MEMORY[DSD_Complex.NAMES[c]], products))
            except KeyError as err:
                raise PilFormatError("Cannot find complex: {}.".format(err))
            anon = Reaction(reactants, products, rtype, rate)

        if anon.rateunits != units:
            raise SystemExit("Rate units must be given in {}, not: {}.".format(
                        reaction.rateunits, units))
        return anon

    raise PilFormatError('unknown keyword: {}'.format(line[0]))

