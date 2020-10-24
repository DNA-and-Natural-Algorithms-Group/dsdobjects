import logging
log = logging.getLogger(__name__)

from .singleton import Singleton, SingletonError
from .utils import flint, convert_units
from .complex_utils import (SecondaryStructureError,
                            make_pair_table, 
                            make_strand_table,
                            strand_table_to_sequence,
                            pair_table_to_dot_bracket,
                            make_loop_index, 
                            wrap,
                            split_complex_pt,
                            rotate_complex_db)

class ObjectInitError(Exception):
    pass

class DomainS(metaclass = Singleton):
    """ Domain object (Singleton).
    
    Each name for the domain will instantiate this class only once.
    Initialization requires a specification of both name and length. Once 
    the domain is initialized [a = Domain('a', 15)], one can access the
    corresponding object or its complement omitting the length 
    [a = Domain('a')]. If, however, a wrong length is provided, then 
    this will raise a DSDObjectsError [b = Domain('a*', 10)]
    """

    @classmethod
    def identifiers(cls, name, length = None):
        """ tuple: A method that must be accessible without initializing the object. """
        return ((name, length), name, {}) if length is not None else (None, name, {})

    def __init__(self, name, length):
        self._name = name
        self._length = length

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
    then this returns the existing complex.

    Sequence and structure can be specified on the domain or on the nucleotide
    level, but they have to be of the same length. Although not implemented,
    one may define special characters for secondary structures that are more
    diverse than a regular dot-bracket string, e.g. 'h' = '(', '.', ')'

    Args:
        sequence (list): A domain-level or nucleotide-level sequence.
        structure (list): A structure in dot-bracket notation.
        name (str, optional): Name of this domain.
    """

    @classmethod
    def identifiers(cls, sequence, structure, name = None, **kwargs):
        """ tuple: A method that must be accessible without initializing the object. """
        if sequence is None or structure is None:
            if name is None:
                raise ObjectInitError('Insufficient arguments for Complex initialization.')
            return (None, name, {})

        if len(sequence) != len(structure):
            raise ObjectInitError(f'Complex initialization: {len(sequence)} != {len(structure)}.')

        cdict = {}
        for e, (rseq, rstr) in enumerate(rotate_complex_db(sequence, structure)):
            rcplx = tuple((tuple(map(str, rseq)), tuple(rstr)))
            if rcplx in cls._instanceCanon:
                canon = rcplx
                turns = e
                break
            cdict[rcplx] = e # How many rotations to the canonical form
        else:
            canon = sorted(cdict, key = lambda x:(x[0], x[1]))[0]
            turns = cdict[canon] # How many rotations to the canonical form

        tot = len(make_strand_table(sequence))
        turns = wrap(-turns, tot) # How many rotations from the canonical form

        newargs = {'canon': canon, 'turns': turns}
        return (canon, name, newargs)

    def __init__(self, sequence, structure, name = None, 
                 canon = None, 
                 turns = None):
        # This must have been set by the identifiers method.
        assert canon is not None
        assert turns is not None

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
        self._enclosed_domains = None
        self._exterior_domains = None
        self._exterior_loops = None

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
        tot = len(self.__strand_table)
        t = wrap(-self._turns + value, tot)
        for e, (seq, sst) in enumerate(self.rotate()):
            if e == t:
                self._sequence = seq
                self._structure = sst
                self._turns = wrap(value, tot)
                break
        else:
            raise ValueError('Something went terribly wrong... ?')

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

    def rotate(self, turns = None):
        return rotate_complex_db(self._sequence, self._structure, turns = turns)
    
    @property
    def sequence(self):
        """ list: sequence the complex object. """
        return self._sequence[:]

    @property
    def strand_table(self):
        return make_strand_table(self._sequence)

    @property
    def structure(self):
        """ list: the complex structure. """
        return self._structure[:]

    @property
    def pair_table(self):
        """ returns a structure in multistranded pair-table format. """
        # Make a new pair_table every time, it might get modified.
        return make_pair_table(self._structure)

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
            self._loop_index, self._exterior_loops = make_loop_index(self.__pair_table)
        return self._loop_index

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
        for si, strand in enumerate(self.__pair_table):
            for di, domain in enumerate(strand):
                loc = (si,di)
                cloc = self._pair_table[si][di] 
                if not (cloc is None or self.get_domain(loc) == ~self.get_domain(cloc)):
                    return False
        return True

    @property
    def is_connected(self):
        if not self._loop_index:
            try :
                self._loop_index, self._exterior_loops = self.__loop_index
            except SecondaryStructureError as e:
                return False
        return True

    def split(self):
        stab = self.__strand_table
        ptab = self.__pair_table
        for st, pt in split_complex_pt(stab, ptab):
            nseq = strand_table_to_sequence(st)
            nsst = pair_table_to_dot_bracket(pt)
            yield self.__class__(nseq, nsst)
        return

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name}, {self.kernel_string})'

    def __str__(self):
        return f'{self.name}'

    def __len__(self):
        raise NotImplementedError('Ambiguous parameter for ComplexS')

    def __eq__(self, other):
        """ Test if two complexes are equal. """
        if not isinstance(other, ComplexS):
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
            assert name in [x.name for x in complexes]
            if name is None:
                name = complexes[0].name
            nargs = {'canon': complexes}
        return (complexes, name, nargs)

    def __init__(self, complexes, name, canon = None):
        self._complexes = complexes
        self._representative = next(x for x in complexes if x.name == name)
        self._canonical_form = canon

    @property
    def complexes(self):
        """ A list of complexes in the resting set. """
        for x in self._complexes:
            yield x

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
        """ Two resting sets are equal if their complexes are equal """
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

    @classmethod
    def identifiers(cls, reactants, products, rtype):
        """ tuple: A method that must be accessible without initializing the object. """
        react = tuple(sorted(list(map(lambda x: x.canonical_form, reactants))))
        prods = tuple(sorted(list(map(lambda x: x.canonical_form, products))))
        canon = tuple((react, prods, rtype))
        name = "reaction [{}] {} -> {}".format(rtype,
                                               " + ".join([x.name for x in reactants]), 
                                               " + ".join([x.name for x in products]))
        newargs = {'canon': canon}
        return (canon, name, newargs)

    def __init__(self, reactants, products, rtype = None, canon = None):
        self._reactants = reactants
        self._products = products
        self._rtype = rtype

        # rate constant in counts per volume
        self._const = None
        self._units = None
        #self._Mol = None
        #self._Mol_unit = None
        #self._counts = None
        #self._count_unit = None
        #self._volume = None
        #self._volume_unit = None
        #self._time = None
        #self._time_unit = None

        assert canon is not None
        self._canonical_form = canon

    @property
    def reactants(self):
        """list: list of reactants. """
        return self._reactants[:]

    @reactants.setter
    def reactants(self, value):
        raise SingletonError(f'{self.__class__.__name__} object reactants is immutable!')

    @property
    def products(self):
        """list: list of products. """
        return self._products[:]

    @products.setter
    def products(self, value):
        raise SingletonError(f'{self.__class__.__name__} object products is immutable!')

    @property
    def rtype(self):
        """str: *peppercorn* reaction type (bind21, condensed, ...) """
        return self._rtype

    @rtype.setter
    def rtype(self, value):
        raise SingletonError(f'{self.__class__.__name__} object rtype is immutable!')

    @property
    def name(self):
        """str: *peppercorn* reaction type (bind21, condensed, ...) """

        if self._const:
            rc = ' = {}{}'.format(flint(self._const), f' {self._units}' if self._units else '')
        else:
            rc = ''
        return "reaction [{}{}] {} -> {}".format(self.rtype, rc,
               " + ".join([x.name for x in self.reactants]), 
               " + ".join([x.name for x in self.products]))

    @property
    def arity(self):
        """(int, int): number of reactants, number of products."""
        return (len(self._reactants), len(self._products))

    @property
    def rate_constant(self):
        return flint(self._const), self._units

    def set_rate_parameters(self, constant, units = None):
        self._const = constant
        self._units = units

    def rateformat(self, output_units):
        """Set reaction rate constant and units."""
        if self._units is None:
            raise ObjectInitError(f'Cannot change the units of the rate constant: {units}.')

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

