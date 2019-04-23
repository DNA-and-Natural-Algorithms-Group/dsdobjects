# python3 
#
# dsdobjects/prototypes.py
#
# Commonly useful DSD object definitions.
#   - functionality from here may be incorporated into base_classes if generally useful.
#   - you should be able to copy that file as is into your project if you need custom changes.
#   - please consider providing thoughts about missing functionality
#

from dsdobjects import clear_memory
from dsdobjects import DL_Domain, SL_Domain, DSD_Complex, DSD_Reaction, DSD_RestingSet, DSD_StrandOrder
from dsdobjects import DSDObjectsError, DSDDuplicationError
from dsdobjects.utils import split_complex

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

    def __init__(self, *kargs, **kwargs):
        super(Domain, self).__init__(*kargs, **kwargs)

    #@property
    #def complement(self):
    #    # If we initialize the complement, we need to know the class.
    #    if self._complement is None:
    #        cname = self._name[:-1] if self.is_complement else self._name + '*'
    #        if cname in DL_Domain.MEMORY:
    #            self._complement = DL_Domain.MEMORY[cname]
    #        else :
    #            self._complement = LogicDomain(cname, self.dtype, self.length)
    #    return self._complement
    
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

    @property
    def pair_table(self):
        return super(Complex, self).pair_table
    
    @pair_table.setter
    def pair_table(self, pt):
        self._pair_table = pt

    def full_string(self):
        return "Complex(%s): %s %s" % (
            self.name, str(self.strands), str(self.structure))

    def get_structure(self, loc):
        return self.get_paired_loc(loc)

    def triple(self, *loc):
        # overwrite standard func
        return (self.get_domain(loc), self.get_paired_loc(loc), loc)

    def get_strand(self, loc):
        """
        Returns the strand at the given index in this complex
        """
        if(loc is not None):
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

    def rotate_location(self, loc, n=None):
        return self.rotate_pairtable_loc(loc, n)

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

    def __cmp__(self, other):
        """
        Two complexes are compared on the basis of their complexes
        """
        return cmp(self.canonical_form, other.canonical_form)

class RestingSet(DSD_RestingSet):
    def __init__(self, *kargs, **kwargs):
        super(RestingSet, self).__init__(*kargs, **kwargs)

    def __str__(self):
        return self.name

    def __len__(self):
        """
        The number of species in a resting set
        """
        return len(self._complexes)

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

        def format_rate_units(rate):
            if time == 's':
                pass
            elif time == 'm':
                rate *= 60
            elif time == 'h':
                rate *= 3600
            else :
                raise NotImplementedError
        
            if molarity == 'M':
                pass
            elif molarity == 'mM':
                if self.arity[0] > 1:
                    factor = self.arity[0] - 1
                    rate /= (factor * 1e3)
            elif molarity == 'uM':
                if self.arity[0] > 1:
                    factor = self.arity[0] - 1
                    rate /= (factor * 1e6)
            elif molarity == 'nM':
                if self.arity[0] > 1:
                    factor = self.arity[0] - 1
                    rate /= (factor * 1e9)
            elif molarity == 'pM':
                if self.arity[0] > 1:
                    factor = self.arity[0] - 1
                    rate /= (factor * 1e12)
            else :
                raise NotImplementedError

            return rate

        rate = format_rate_units(self.rate) if self.rate else float('nan')
        units = "/{}".format(molarity) * (self.arity[0] - 1) + "/{}".format(time)

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

        Returns:
            StrandOrder, pairtable-reactants, pairtable-products.

        """
        # find a common pairtable representation for input and output

        # get the StrandOrder of inputs and the StrandOrder of outputs
        # e.g.:
        #   (sX + sY) + sZ => (sY + sZ) + sX
        # we know that all intermediates must be PK free, does that guarantee that such an order exists?
        #
        #   (X+Y) + Z => (XZY) => (ZY) + X => YZ + X

        # If len(reactants) == 1 or len(products)==1, then we have the strand order.
        # else it can be very complicated..., but frankly, we just have to try...

        so = None  # The common strand order.
        pt1 = None # Pair table of reactants
        pt2 = None # Pair table of products

        rotations = 0
        if len(self.reactants) == 1:
            cplx = self.reactants[0]
            if isinstance(cplx, RestingSet):
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
            if isinstance(cplx, RestingSet):
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

        if any(map(lambda c: isinstance(c, RestingSet), cxs)):
            cxs = list(map(lambda x: x.canonical, cxs))

        assert len(cxs) == 2
        for rot1 in cxs[0].rotate():
            for rot2 in cxs[1].rotate():
                rotations = None
                try:
                    so2 = StrandOrder(rot1.sequence + ['+'] + rot2.sequence)
                    #rotations = so2.rotations
                except DSDDuplicationError as e : 
                    so2 = e.existing
                    rotations = e.rotations
                if so == so2:
                    fake = Complex(rot1.sequence + ['+'] + rot2.sequence, rot1.structure + ['+'] + rot2.structure)
                    if rotations:
                        for x in range(len(so)-rotations):
                            fake = fake.rotate_once()
                    if pt1: pt2 = fake.pair_table
                    else :  pt1 = fake.pair_table

        # What have we got?
        return (so, pt1, pt2)

class StrandOrder(DSD_StrandOrder):
    pass

    #def __new__(cls, sequence, name='', prefix='StrandOrder_', memorycheck=True):
    #    # The new method returns the present instance of an object, if it exists
    #    self = DSD_StrandOrder.__new__(cls)
    #    try:
    #        super(StrandOrder, self).__init__(sequence, name, prefix, memorycheck)
    #    except DSDDuplicationError as e :
    #        return e.existing
    #    return self

    #def __init__(self, sequence, name='', prefix='StrandOrder_', memorycheck=True):
    #    # Remove default initialziation to get __new__ to work
    #    pass

