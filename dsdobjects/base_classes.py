import logging
log = logging.getLogger(__name__)

from itertools import chain

from .singleton import Singleton, SingletonError, show_singletons
from .utils import flint, convert_units
from .complex_utils import (SecondaryStructureError,
                            make_pair_table, 
                            make_strand_table,
                            strand_table_to_sequence,
                            pair_table_to_dot_bracket,
                            make_loop_index, 
                            wrap,
                            split_complex_pt,
                            rotate_complex_once)

class ObjectInitError(Exception):
    pass

def show_memory():
    """ Shows all objects of base_classes.py that are in use.

    This function is for debugging only! The memory of a Singleton clears
    automatically if there is no hard reference to the object. However, one has
    to make -- possibly unexpected, un-pythonic -- handstands to ensure that
    hard references are cleaned up. For example when working with dictionaries
    and sets, where the "del" is not enough, but you need to clear() all
    entries, before you return the results of a function.
    """
    for x in chain(show_singletons(DomainS),
                   show_singletons(StrandS),
                   show_singletons(ComplexS),
                   show_singletons(MacrostateS),
                   show_singletons(ReactionS)):
        yield x

class DomainS(metaclass = Singleton):
    """ Domain object (Singleton).
    
    Each name for the domain will instantiate this class only once.
    Initialization requires a specification of both name and length. Once 
    the domain is initialized [a = Domain('a', 15)], one can access the
    corresponding object or its complement omitting the length 
    [a = Domain('a')]. If, however, a wrong length is provided, then 
    this will raise a DSDObjectsError [b = Domain('a*', 10)]
    """
    DTYPE_CUTOFF = 8 
    SHORT_DOM_LEN = 5
    LONG_DOM_LEN = 15
    PREFIX = 'd'
    ID = 1

    @classmethod
    def identifiers(cls, name = None, length = None, prefix = None, dtype = None):
        """ tuple: A method that must be accessible without initializing the object. """
        if name is None:
            if prefix is None:
                prefix = f'{cls.PREFIX}'
            name = f'{prefix}{cls.ID}'
        if length is None:
            length = cls.SHORT_DOM_LEN if dtype == 'short' else \
                     cls.LONG_DOM_LEN if dtype == 'long' else None
        elif dtype and not (dtype == 'short') == (length <= cls.DTYPE_CUTOFF):
            raise ObjectInitError(f'Conflicting arguments {dtype} and {length}.')

        cname = f'{name[:-1]}' if name[-1] == '*' else f'{name}*' 
        newargs = {}
        if length is None and name[-1] == '*':
            # Allow initialization of a complementary domain without length!
            try:
                length = len(cls(cname, length = None))
                newargs = {'length': length}
            except SingletonError:
                pass
        elif length and name[-1] != '*':
            # Forbid initialization of a non-complementary domain with conflicting length.
            clength = length
            try:
                clength = len(cls(cname))
                cls(cname, length = length)
            except SingletonError:
                if clength != length:
                    raise SingletonError(f'Duplicate Singleton {cls.__name__}: name ({name}) has the wrong length {length} vs {clength}!')
        return ((name, length), name, newargs) if length is not None else (None, name, {})

    def __init__(self, name = None, length = None, prefix = None, dtype = None):
        if name is None:
            if prefix is None:
                prefix = f'{self.__class__.PREFIX}'
            name = f'{prefix}{self.__class__.ID}'
            self.__class__.ID += 1
        if length is None:
            length = self.__class__.SHORT_DOM_LEN if dtype == 'short' else \
                     self.__class__.LONG_DOM_LEN if dtype == 'long' else None
        self._name = name
        self._length = length
        self.sequence = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        raise SingletonError(f'{self.__class__.__name__} object name is immutable!')

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        raise SingletonError(f'{self.__class__.__name__} object length is immutable!')

    @property
    def dtype(self):
        return 'short' if self.length <= self.__class__.DTYPE_CUTOFF else 'long'

    @property
    def canonical_form(self):
        return self.__class__.identifiers(self.name, self.length)

    @property
    def is_complement(self):
        return self.name[-1] == '*'

    @property
    def cname(self):
        """ str: the name of the complementary domain. """
        return self.name[:-1] if self.is_complement else self.name + '*'

    @property
    def complement(self):
        """ obj: the complementary domain object. """
        return self.__class__(self.cname, self.length)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name}, {self.length})"

    def __str__(self):
        return f'{self.name}'

    def __len__(self):
        return self.length

    def __invert__(self): 
        return self.complement

    def __eq__(self, other):
        """ Test if two domains are equal. """
        # We use DomainS here, not self.__class__!
        # Equality is based on the canonical form, otherwise use "is"
        if not isinstance(other, DomainS):
            return False
        return (self.name, self.length) == (other.name, other.length)

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    def __le__(self, other):
        return self.name <= other.name

    def __ge__(self, other):
        return self.name >= other.name

    def __hash__(self):
        return hash(self.name)

class ComplexS(metaclass = Singleton):
    """ Complex object (Singleton).

    If the same complex is initialized twice (e.g. in a different rotation),
    then this returns the existing complex.  Sequence and structure should be
    specified on the domain level, and they have to be of the same length. 

    Args:
        sequence (list): A domain-level sequence.
        structure (list): A structure in dot-bracket notation.
        name (str, optional): A name tag unique for this complex.
    """
    PREFIX = 'c'
    ID = 1

    @classmethod
    def identifiers(cls, sequence, structure, name = None, prefix = None, **kwargs):
        """ tuple: A method that must be accessible without initializing the object. """
        if sequence is None:
            if name is None:
                raise ObjectInitError('Insufficient arguments for Complex initialization.')
            canon = None
            newargs = {}
        else:
            if name is None:
                name = f'{cls.PREFIX}{cls.ID}' if prefix is None else f'{prefix}{cls.ID}'
            if len(sequence) != len(structure):
                raise ObjectInitError('Complex initialization error: ' + \
                                     f'{len(sequence)} != {len(structure)}.')
            cdict = {} # Find canonical form.
            rseq, rstr = sequence, structure
            for e in range(len(make_strand_table(sequence))):
                rcplx = tuple((tuple(map(str, rseq)), tuple(rstr)))
                if rcplx in cls._instanceCanon:
                    canon = rcplx
                    turns = e
                    break
                cdict[rcplx] = e # How many rotations to the canonical form
                rseq, rstr = rotate_complex_once(rseq, rstr)
            else:
                canon = sorted(cdict, key = lambda x:(x[0], x[1]))[0]
                turns = cdict[canon] # How many rotations to the canonical form
            tot = len(make_strand_table(sequence))
            turns = wrap(-turns, tot) # How many rotations from the canonical form
            newargs = {'canon': canon, 'turns': turns, 'rcplxs': cdict.keys()}
        return (canon, name, newargs)

    def __init__(self, sequence, structure, name = None, 
                 prefix = None, canon = None, turns = None, rcplxs = None):
        # This must have been set by the identifiers method.
        cls = self.__class__
        assert canon is not None
        assert turns is not None
        if name is None:
            name = f'{cls.PREFIX}{cls.ID}' if prefix is None else f'{prefix}{cls.ID}'
            self.__class__.ID += 1

        # Private variables:
        self._sequence = sequence
        self._structure = structure
        self._name = name
        self._canon = canon
        self._turns = turns

        # Initialized on demand:
        self._strand_table = None
        self._pair_table = None
        self._loop_index = None
        self._domains = None
        self._exterior_domains = None
        self._enclosed_domains = None
        self._exterior_loops = None
        self._concentration = None

        # A speedup that uses some memory ...
        for rcplx in rcplxs:
            cls._instanceCanon[rcplx] = self

    @property
    def name(self):
        """ str: name of the complex object. """
        return self._name

    @name.setter
    def name(self, value):
        raise SingletonError(f'{self.__class__.__name__} object name is immutable!')

    @property
    def canonical_form(self):
        """ tuple: lexicographically unique sorting of sequence & structure. """
        return self._canon

    @canonical_form.setter
    def canonical_form(self, value):
        raise SingletonError(f'{self.__class__.__name__} object canonical_form is immutable!')

    @property
    def turns(self):
        """ Number of cyclic permutations from canonical form to representation. """
        return self._turns

    @turns.setter
    def turns(self, value):
        # Turns = 0 rotates the object into the canonical form.
        # Turns = 1 rotates the object into canonical form + 1 turn.
        tot = self.size
        t = wrap(-self._turns + value, tot)
        for e, (seq, sst) in enumerate(self.rotate()):
            if e == t:
                self._sequence = seq
                self._structure = sst
                self._turns = wrap(value, tot)
                break
        else:
            raise ObjectInitError('Something went terribly wrong when rotating the complex.')

    @property
    def sequence(self):
        """ list: sequence the complex object. """
        return iter(self._sequence)

    @property
    def strand_table(self):
        if not self._strand_table:
            self._strand_table = make_strand_table(self._sequence)
        for strand in self._strand_table:
            yield list(strand) # Deepcopy

    @property
    def structure(self):
        """ list: the complex structure. """
        return iter(self._structure)

    @property
    def pair_table(self):
        """ returns a structure in multistranded pair-table format. """
        if not self._pair_table:
            self._pair_table = make_pair_table(self._structure)
        for locs in self._pair_table:
            yield list(locs) # Deepcopy

    @property
    def __strand_table(self):
        if not self._strand_table:
            self._strand_table = make_strand_table(self._sequence)
        return self._strand_table

    @property
    def __pair_table(self):
        if not self._pair_table:
            self._pair_table = make_pair_table(self._structure)
        return self._pair_table

    @property
    def __loop_index(self):
        if not self._loop_index:
            self._loop_index, self._exterior_loops = make_loop_index(self.pair_table)
        return self._loop_index

    @property
    def size(self):
        return len(self.__strand_table)

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
            self._concentration = (mode, value, unit)

    def concentrationformat(self, out):
        mod = self._concentration[0]
        val = self._concentration[1]
        uni = self._concentration[2]
        val = convert_units(val, uni, out)
        return (mod, val, out)

    def rotate(self, turns = None):
        """ Returns every rotation of the sequence, structure pair for the complex.

        It starts with the default representation.
        """
        if turns is None:
            turns = self.size
        yield list(self.sequence), list(self.structure)
        x, y = self._sequence, self._structure
        for _ in range(turns-1):
            x, y = rotate_complex_once(x, y)
            yield x, y

    def rotate_pt(self, turns = None):
        """ A wrapper for rotate() which returns strand table and pair table. """
        for (x, y) in self.rotate(turns):
            yield make_strand_table(x), make_pair_table(y)
    
    # ------ can be mutable but must yield the same canonical form!
    def strand_length(self, pos):
        return len(self.__strand_table[pos])
 
    def get_loop_index(self, loc):
        return self.__loop_index[loc[0]][loc[1]]

    def get_domain(self, loc):
        return self.__strand_table[loc[0]][loc[1]]
    
    def get_paired_loc(self, loc):
        """ 
        Returns the paired element in the pair-table. 
        Raises: IndexError if there are negative elements in loc
        """
        if loc[0] < 0 or loc[1] < 0:
            raise IndexError
        return self.__pair_table[loc[0]][loc[1]]

    def rotate_pairtable_loc(self, loc, n):
        """ Maps the locus of a given pair-table to a new rotation.  """
        return (wrap(loc[0] - n, self.size), loc[1])

    @property
    def domains(self):
        if not self._domains:
            self._domains = set(self.sequence)
            self._domains -= set('+')
        return self._domains

    @property
    def enclosed_domains(self):
        if not self._enclosed_domains:
            _ = self.exterior_domains
        return self._enclosed_domains

    @property
    def exterior_domains(self):
        """
        Returns all domains in exterior loops.
        """
        if not self._exterior_domains:
            self._exterior_domains = []
            self._enclosed_domains = []
            for si, strand in enumerate(self.__loop_index):
                for di, domain in enumerate(strand):
                    if self._loop_index[si][di] in self._exterior_loops:
                        if self._pair_table[si][di] is None:
                            self._exterior_domains.append((si, di))
                    elif self._pair_table[si][di] is None:
                            self._enclosed_domains.append((si, di))
        return self._exterior_domains

    # Sanity Checks
    @property
    def is_domainlevel_complement(self):
        """
        Determines whether the structure includes pairs only between complementary domains.
        Returns True if all paired domains are complementary, raises an Exception otherwise
        """
        for si, strand in enumerate(self.pair_table):
            for di, domain in enumerate(strand):
                loc = (si,di)
                cloc = self._pair_table[si][di] 
                if not (cloc is None or self.get_domain(loc) == ~self.get_domain(cloc)):
                    return False
        return True

    @property
    def is_connected(self):
        if not self._loop_index:
            try:
                _ = self.__loop_index
            except SecondaryStructureError as e:
                return False
        return True

    def split(self):
        stab = self.__strand_table
        ptab = self.__pair_table
        for st, pt in split_complex_pt(stab, ptab):
            nseq = strand_table_to_sequence(st)
            nsst = pair_table_to_dot_bracket(pt)
            try:
                yield self.__class__(nseq, nsst)
            except SingletonError as err:
                if err.existing is None:
                    log.warning(f'Automated naming of {self.__class__} object failed. You may have to change {self.__class__}.PREFIX.')
                    raise err
                yield err.existing
        return

    @property
    def kernel_string(self):
        """ str: print sequence and structure in `kernel` notation. """
        seq = self._sequence
        sst = self._structure
        knl = ''
        for i in range(len(seq)):
            if sst[i] == '+':
                knl += str(sst[i]) + ' '
            elif sst[i] == ')':
                knl += str(sst[i]) + ' '
            elif sst[i] == '(':
                knl += str(seq[i]) + str(sst[i]) + ' '
            else:
                knl += str(seq[i]) + ' '
        return knl[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name}, {self.kernel_string})'

    def __str__(self):
        return f'{self.name}'

    def __len__(self):
        raise NotImplementedError('Ambiguous parameter for ComplexS')

    def __eq__(self, other):
        """ Test if two complexes are equal. """
        # We use ComplexS here, not self.__class__!
        # Equality is based on the canonical form, otherwise use "is"
        if not isinstance(other, ComplexS):
            return False
        return self.canonical_form == other.canonical_form

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        assert isinstance(other, ComplexS)
        return self.canonical_form < other.canonical_form

    def __gt__(self, other):
        assert isinstance(other, ComplexS)
        return self.canonical_form > other.canonical_form

    def __le__(self, other):
        assert isinstance(other, ComplexS)
        return self.canonical_form <= other.canonical_form

    def __ge__(self, other):
        assert isinstance(other, ComplexS)
        return self.canonical_form >= other.canonical_form

    def __hash__(self):
        return hash(self.canonical_form)

class StrandS(ComplexS):
    PREFIX = 's'
    ID = 1
    @classmethod
    def identifiers(cls, sequence, name = None, prefix = None, **kwargs):
        """ tuple: A method that must be accessible without initializing the object. """
        if sequence is None:
            if name is None:
                raise ObjectInitError('Insufficient arguments for Strand initialization.')
            canon = None
            newargs = {}
        elif '+' in sequence:
            raise NotImplementedError('ComplexS "strand" mode must only contain a single strand.')
        else:
            if name is None:
                name = f'{cls.PREFIX}{cls.ID}' if prefix is None else f'{prefix}{cls.ID}'
            sstr = tuple('*' for _ in range(len(sequence)))
            canon = (tuple(map(str, sequence)), sstr)
            newargs = {'canon': canon, 'turns': 0}
        return (canon, name, newargs)

    def __init__(self, sequence, name = None, prefix = None, canon = None, turns = None):
        # This must have been set by the identifiers method.
        cls = self.__class__
        assert turns == 0
        assert canon is not None
        if name is None:
            name = f'{cls.PREFIX}{cls.ID}' if prefix is None else f'{prefix}{cls.ID}'
            self.__class__.ID += 1

        # Private variables:
        self._sequence = sequence
        self._name = name
        self._structure = None
        self._canon = canon
        self._turns = turns

        # Initialized on demand:
        self._strand_table = None
        self._domains = None
        self._concentration = None

    @property
    def structure(self):
        return None

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name}, {" ".join(map(str, self._sequence))})'

class MacrostateS(metaclass = Singleton):
    """ A set of complexes (singleton). 

    Macrostates are initialized with a name, where the name points to a
    particular complex. 
    """
    @classmethod
    def identifiers(cls, complexes = None, name = None):
        """ tuple: A method that must be accessible without initializing the object. """
        if complexes is None:
            assert name is not None
            nargs = {}
        else:
            complexes = tuple(sorted(complexes, key = lambda x: x.canonical_form))
            nargs = {'canon': complexes}
            if name is None:
                name = complexes[0].name
                nargs['name'] = name
            else:
                assert name in [x.name for x in complexes]
        return (complexes, name, nargs)

    def __init__(self, complexes, name, canon = None):
        self._complexes = complexes
        self._representative = next(x for x in complexes if x.name == name)
        self._canonical_form = canon

    @property
    def complexes(self):
        """ A list of complexes in the resting set. """
        return iter(self._complexes)

    @complexes.setter
    def complexes(self, value):
        raise SingletonError(f'{self.__class__.__name__} object complexes is immutable!')

    @property
    def representative(self):
        return self._representative

    @representative.setter
    def representative(self, value):
        raise SingletonError(f'{self.__class__.__name__} object repersentative is immutable!')

    @property
    def canonical_form(self):
        return self._canonical_form

    @property
    def name(self):
        return self.representative.name

    @property
    def kernel_string(self):
        return self.representative.kernel_string

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name}, {[x.name for x in self.complexes]})"

    def __str__(self):
        return f'{self.name}'

    def __len__(self):
        """ The number of species in a macrostate """
        return len(self._complexes)

    def __eq__(self, other):
        # We use MacrostateS here, not self.__class__!
        # Equality is based on the canonical form, otherwise use "is"
        if not isinstance(other, MacrostateS):
            return False
        return (self.canonical_form == other.canonical_form)

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        """ ReactionPathway objects are sorted by their canonical form.  """
        return (self.canonical_form < other.canonical_form)

    def __gt__(self, other):
        """ ReactionPathway objects are sorted by their canonical form.  """
        return (self.canonical_form > other.canonical_form)

    def __le__(self, other):
        return self.canonical_form <= other.canonical_form

    def __ge__(self, other):
        return self.canonical_form >= other.canonical_form

    def __hash__(self):
        return hash(self.canonical_form)

class ReactionS(metaclass = Singleton):
    """ A reaction between complexes or macrostates (singleton). 

    Args:
      reactants (list): A list of reactants. Reactants can be 
        :obj:`DSD_Macrostate()` or :obj:`DSD_Complex()` objects.
      products (list): A list of products. Products can be
        :obj:`DSD_Macrostate()` or :obj:`DSD_Complex()` objects.
      rtype (str, optional): Reaction type, e.g. bind21, condensed, .. Defaults to None.
      rate (flt, optional): Reaction rate. A reaction rate 
    """
    RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    @classmethod
    def identifiers(cls, reactants, products, rtype, name = None):
        """ tuple: A method that must be accessible without initializing the object. """
        if name is not None and reactants is None and products is None and rtype is None:
            return (None, name, {})
        react = tuple(sorted([x.canonical_form for x in reactants]))
        prods = tuple(sorted([x.canonical_form for x in products]))
        canon = tuple((react, prods, rtype))
        newargs = {'canon': canon}
        if name is None:
            name = "[{}] {} -> {}".format(rtype,
                    " + ".join([x.name for x in sorted(reactants, 
                                                       key = lambda y: y.canonical_form)]), 
                    " + ".join([x.name for x in sorted(products, 
                                                       key = lambda y: y.canonical_form)]))
            newargs['name'] = name
        return (canon, name, newargs)

    def __init__(self, reactants, products, rtype, name = None, canon = None):
        self._reactants = sorted(reactants, key = lambda y: y.canonical_form)
        self._products = sorted(products, key = lambda y: y.canonical_form)
        self._rtype = rtype

        # rate constant in counts per volume
        self._const = None
        self._units = None

        assert name is not None
        self._name = name
        assert canon is not None
        self._canonical_form = canon

    @property
    def reactants(self):
        """list: list of reactants. """
        return iter(self._reactants)

    @reactants.setter
    def reactants(self, value):
        raise SingletonError(f'{self.__class__.__name__} object reactants is immutable!')

    @property
    def products(self):
        """list: list of products. """
        return iter(self._products)

    @products.setter
    def products(self, value):
        raise SingletonError(f'{self.__class__.__name__} object products is immutable!')

    @property
    def rtype(self):
        """str: reaction type (bind21, condensed, ...) """
        return self._rtype

    @rtype.setter
    def rtype(self, value):
        raise SingletonError(f'{self.__class__.__name__} object rtype is immutable!')

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        raise SingletonError(f'{self.__class__.__name__} object rtype is immutable!')

    @property
    def reaction_string(self):
        if self._const:
            rc = ' = {:10g}{:5s} '.format(flint(self._const), f' {self._units}' if self._units else '')
        else:
            rc = ''
        return "reaction [{:12s}{:12s}] {} -> {}".format(
                self.rtype, rc, " + ".join([x.name for x in self.reactants]), 
                                " + ".join([x.name for x in self.products]))
    @property
    def arity(self):
        """ (int, int): number of reactants, number of products."""
        return len(self._reactants), len(self._products)

    @property
    def rate_constant(self):
        if self._const is None:
            return None, None
        else:
            return flint(self._const), self._units

    @rate_constant.setter
    def rate_constant(self, tup):
        if isinstance(tup, tuple):
            assert 1 <= len(tup) <= 2
            (constant, units) = tup if len(tup) == 2 else (tup, None)
        else:
            (constant, units) = (tup, None)
        self._const = constant
        self._units = units

    def rateformat(self, output_units):
        """ Set reaction rate constant and units. """
        if self._units is None:
            raise ObjectInitError(f'Cannot change the units of the rate constant: {self._units}.')
        old = self._units.split('/')[1:]
        if len(old) != len(self._reactants):
            raise NotImplementedError(f'Cannot interpret the format of units: {self._units}')
        new = output_units.split('/')[1:]
        if len(new) != len(self._reactants):
            raise NotImplementedError(f'Cannot interpret the format of units: {output_units}')

        newc = self._const
        for i, o in zip(old, new):
            newc = convert_units(newc, o, i) # 1/M 1/s
        return newc, output_units

    @property
    def kernel_string(self):
        return "{} -> {}".format("  +  ".join(r.kernel_string for r in self.reactants),
                                 "  +  ".join(p.kernel_string for p in self.products))

    @property
    def canonical_form(self):
        return self._canonical_form

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name})"

    def __str__(self):
        return f'{self.name}'

    def __eq__(self, other):
        # We use ReactionS here, not self.__class__!
        # Equality is based on the canonical form, otherwise use "is"
        if not isinstance(other, ReactionS):
            return False
        return self.canonical_form == other.canonical_form

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.canonical_form < other.canonical_form

    def __gt__(self, other):
        return self.canonical_form > other.canonical_form

    def __le__(self, other):
        return self.canonical_form <= other.canonical_form

    def __ge__(self, other):
        return self.canonical_form >= other.canonical_form

    def __hash__(self):
        return hash(self.canonical_form)

