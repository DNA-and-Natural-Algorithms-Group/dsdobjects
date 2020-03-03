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

import logging
log = logging.getLogger(__name__)

from collections import namedtuple

from dsdobjects.core import DSDObjectsError, DSDDuplicationError
from dsdobjects.core import SequenceConstraint # just pass it on ...
from dsdobjects.core import DL_Domain, SL_Domain 
from dsdobjects.core import DSD_Complex, DSD_Reaction, DSD_Macrostate, DSD_StrandOrder
from dsdobjects.utils import split_complex, natural_sort, convert_units

class LogicDomain(DL_Domain):
    """
    Represents a single domain. We allow several options for specifying domain
    properties. Domains might have an explicit integer (bp) length, or may be
    designated as short or long. If the latter method is used, the code will use
    the relevant constant as the integer domain length.
    """
    
    def __init__(self, name='', prefix='d', dtype=None, length=None):
        # Assign name
        if name == '':
            if prefix == '':
                raise NuskellObjectError('NuskellDomain prefix must not be empty!')
            elif prefix[-1].isdigit():
                raise NuskellObjectError('NuskellDomain must not end with a digit!')
            name = prefix + str(NuskellDomain.ID)
            NuskellDomain.ID += 1
        super(NuskellDomain, self).__init__(name, dtype, length)



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

        self.nucleotides = None
        return self

    def __init__(self, name, dtype = None, length = None):
        # Remove default initialziation to get __new__ to work
        pass

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

    CONCENTRATION = namedtuple('concentration', 'mode value unit')

    @staticmethod
    def clear_memory(memory=True, names=True, ids=True):
        if memory:
            DSD_Complex.MEMORY = dict()
        if names:
            DSD_Complex.NAMES = dict()
        if ids:
            DSD_Complex.ID = dict()

    def __init__(self, sequence, structure, name='', prefix='cplx', memorycheck=True):
        super(Complex, self).__init__(sequence, structure, name, prefix, memorycheck)
        self._concentration = None # e.g. (initial, 5, nM)
        assert self.is_domainlevel_complement

    @property
    def concentration(self):
        return self._concentration

    @concentration.setter
    def concentration(self, trip):
        if trip is None:
            self._concentration = None
        else:
            (mode, value, unit) = trip
            assert isinstance(value, (int, float))
            mode = 'initial' if mode[0] == 'i' else 'constant'
            self._concentration = Complex.CONCENTRATION(mode, value, unit)

    def concentrationformat(self, out):
        # const = Complex.concentrationformat('M').value
        mod = self._concentration.mode
        val = self._concentration.value
        uni = self._concentration.unit
        val = convert_units(val, uni, out) 
        return Complex.CONCENTRATION(mod, val, out)

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
                    " + ".join(natural_sort(map(str, self.reactants))), " + ".join(natural_sort(map(str, self.products))))
        else :
            return '[{:12g} {:4s} ] {} -> {}'.format(rate, units,
                    " + ".join(natural_sort(map(str, self.reactants))), " + ".join(natural_sort(map(str, self.products))))

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
                cplx = cplx.canonical_complex
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
                cplx = cplx.canonical_complex
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
            cxs = list(map(lambda x: x.canonical_complex, cxs))

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


