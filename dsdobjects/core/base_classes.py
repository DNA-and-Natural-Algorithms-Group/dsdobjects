#
# dsdobjects/core/base_classes.py
#   - copy and/or modify together with tests/test_base_classes.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
#  Contributions:
#  This file contains adapted code from various related Python packages
#  coded in the "DNA and Natural Algorithms Group", Caltech:
#   * "DNAObjecs" coded by Joseph Berleant and Joseph Schaeffer 
#   * "Peppercornenumerator" coded by Kathrik Sarma, Casey Grun and Erik Winfree
#   * "Nuskell" coded by Seung Woo Shin
#
# Distributed under the MIT License, use at your own risk.
#

# Python 3 compatibility
from __future__ import absolute_import, division, print_function

from collections import namedtuple

import dsdobjects.utils as utils

class DSDObjectsError(Exception):
    """
    dnaobjects.base_classes error class.
    """

    def __init__(self, message, *kargs):
        if kargs:
            self.message = "{} [{}]".format(message, ', '.join(map(str,kargs)))
        else :
            self.message = message
        super(DSDObjectsError, self).__init__(self.message)

    @property
    def existing(self):
        """return the existing object."""
        return self._existing

    @existing.setter
    def existing(self, value):
        """specify the existing object."""
        self._existing = value

class DSDDuplicationError(Exception):
    """
    dnaobjects.base_classes duplication error class. 
    
    Unique class to return an existing equivalent object from memory, which may
    be used directly by child classes, instead of initializing a new object.
    """
    def __init__(self, message, *kargs):
        if kargs:
            self.message = "{} [{}]".format(message, ', '.join(map(str,kargs)))
        else :
            self.message = message
        self._existing = None
        self._rotations = None
        super(DSDDuplicationError, self).__init__(self.message)

    @property
    def existing(self):
        """return the existing object."""
        return self._existing

    @existing.setter
    def existing(self, value):
        """specify the existing object."""
        self._existing = value

    @property
    def rotations(self):
        """
        :int: return number of rotations to translate a duplicated complex into the
        existing one.  
        """
        return self._rotations

    @rotations.setter
    def rotations(self, value):
        """Set number of rotations to translate a duplicated complex into the
        existing one.

        Args:
            value(int): Number of rotations.
        """
        self._rotations = value

class SequenceConstraint(object):
    """ 
    A nucleic acid sequence constraint in IUPAC format.

    Args:
        sequence(list or str):  A sequence of nucleotides.
        molecule(string, optional): Choice of "DNA" or "RNA" alphabet. Defaults
            to "DNA".
    """
    def __init__(self, sequence, molecule='DNA'):
        assert molecule == 'DNA' or molecule == 'RNA'
        self.ToU = 'T' if molecule == 'DNA' else 'U'

        self._molecule = molecule
        self._sequence = list(sequence)

    @property
    def constraint(self):
        """:str: constraint."""
        return ''.join(self._sequence)

    @property
    def complement(self):
        """:str: complement including wobble base pairs."""
        return ''.join(map(self._iupac_complement, self._sequence))

    @property
    def wc_complement(self):
        """:str: Watson-Crick complement."""
        return ''.join(map(self._wc_complement, self._sequence))

    @property
    def reverse_complement(self):
        """:str: reverse complement including wobble base pairs."""
        return ''.join(map(self._iupac_complement, self._sequence[::-1]))

    @property
    def reverse_wc_complement(self):
        """:str: reverse Watson-Crick complement."""
        return ''.join(map(self._wc_complement, self._sequence[::-1]))

    def add_constraint(self, con):
        """Apply nucleotide sequence constraint and check if it is compatible.

        Args:
          con (str): A new sequence constraint.

        Raises:
          DSDObjectsError: Sequence constraints differ in length.
          DSDObjectsError: Sequence constraints cannot be satisfied.
        """
        con = list(con) # just to make sure...

        if len(self._sequence) != len(con):
            raise DSDObjectsError("Sequence constraints differ in length.")

        new = self._merge_constraints(self._sequence, con)

        if '' in new:
            raise DSDObjectsError("Sequence constraints cannot be satisfied.")
        else:
            self._sequence = new

    def __add__(self, other):
        assert isinstance(other, SequenceConstraint)
        assert self._molecule == other._molecule
        new = SequenceConstraint(self.constraint, self._molecule)
        new.add_constraint(other.constraint)
        return new

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __eq__(self, other):
        return self.constraint == other.constraint

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return str(self.constraint)

    def __len__(self):
        return len(self._sequence)

    def __invert__(self):
        return SequenceConstraint(self.reverse_complement)

    def _merge_constraints(self, con, con2):
        """Return a new list of unified constraints. """
        return list(map(self._iupac_union, zip(con, con2)))

    def _iupac_to_bases(self, nuc):
        # This is more of educational purpose, use _iupac_bin and _bin_iupac
        T = self.ToU
        iupac_dict = {
            'A': 'A',
            'C': 'C',
            'G': 'G',
             T :  T,
            'R': 'AG',      # purine
            'Y': 'C' + T,   # pyrimidine
            'S': 'CG',      # strong
            'M': 'AC',
            'W': 'A' + T,   # weak
            'K': 'G' + T,
            'V': 'ACG',     # not T
            'H': 'AC' + T,  # not G
            'D': 'AG' + T,  # not C
            'B': 'CG' + T,  # not A
            'N': 'ACG' + T}
        return iupac_dict[nuc]

    def _iupac_complement(self, nuc):
        T = self.ToU
        neighbor_dict = {   # ACGT => ACGT
            'A': T,         # 1000 => 0001
            'C': 'G',       # 0100 => 0010
            'G': 'Y',       # 0010 => 0101
             T : 'R',       # 0001 => 1010
            'T': 'R',       # 0001 => 1010
            'R': 'Y',       # 1010 => 0101
            'Y': 'R',       # 0101 => 1010
            'S': 'B',       # 0110 => 0111
            'M': 'K',       # 1100 => 0011
            'W': 'D',       # 1001 => 1011
            'K': 'N',       # 0011 => 1111
            'V': 'B',       # 1110 => 0111
            'H': 'D',       # 1101 => 1011
            'D': 'N',       # 1011 => 1111
            'B': 'N',       # 0111 => 1111
            'N': 'N'}       # 1111 => 1111
        return neighbor_dict[nuc]

    def _wc_complement(self, nuc):
        T = self.ToU
        neighbor_dict = {   # ACGT => ACGT
            'A': T,         # 1000 => 0001
            'C': 'G',       # 0100 => 0010
            'G': 'C',       # 0010 => 0100
             T : 'A',       # 0001 => 1000
            'N': 'N'}       # 1111 => 1111
        return neighbor_dict[nuc]

    def _iupac_union(self, nucs):
        # Return the maximal common constraint
        u = 'N'
        for n in nucs:
            u = self._bin_iupac(self._iupac_bin(u) & self._iupac_bin(n))
        return u

    def _iupac_bin(self, nuc):
        T = self.ToU
        iupac_bin_dict = {  # ACGT
            'A': 8,         # 1000,
            'C': 4,         # 0100,
            'G': 2,         # 0010,
             T : 1,         # 0001,
            'R': 10,        # 1010,  # purine
            'Y': 5,         # 0101,  # pyrimidine
            'S': 6,         # 0110,
            'M': 12,        # 1100,
            'W': 9,         # 1001,
            'K': 3,         # 0011,
            'V': 14,        # 1110,  # not T
            'H': 13,        # 1101,  # not G
            'D': 11,        # 1011,  # not C
            'B': 7,         # 0111,  # not A
            'N': 15}        # 1111,
        return iupac_bin_dict[nuc]

    def _bin_iupac(self, nuc):
        T = self.ToU
        bin_iupac_dict = [  # ACGT
            '',             # 0000  0
            T,              # 0001  1
            'G',            # 0010  2
            'K',            # 0011  3
            'C',            # 0100  4
            'Y',            # 0101  5
            'S',            # 0110  6
            'B',            # 0111  7
            'A',            # 1000  8
            'W',            # 1001  9
            'R',            # 1010 10
            'D',            # 1011 11
            'M',            # 1100 12
            'H',            # 1101 13
            'V',            # 1110 14
            'N']            # 1111 15
        return bin_iupac_dict[nuc]

def clear_memory():
    """ 
    Reset all dsdobjects base-class variables.
    """
    DL_Domain.MEMORY = dict()
    SL_Domain.MEMORY = dict()
    DSD_StrandOrder.MEMORY = dict()
    DSD_StrandOrder.NAMES = dict()
    DSD_StrandOrder.ID = 0
    DSD_Complex.MEMORY = dict()
    DSD_Complex.NAMES = dict()
    DSD_Complex.ID = 0
    DSD_Macrostate.NAMES = dict()
    DSD_Macrostate.MEMORY = dict()
    DSD_Reaction.MEMORY = dict()

class ABC_Domain(object):
    """Abstract base class - domain"""

    def __init__(self):
        super(ABC_Domain, self).__init__()

    @property
    def name(self):
        raise NotImplementedError

    @property
    def complement(self):
        raise NotImplementedError

    @property
    def length(self):
        raise NotImplementedError

    @property
    def sequence(self):
        raise NotImplementedError

    @property
    def constraint(self):
        raise NotImplementedError

    def __len__(self):
        return self.length

    def __eq__(self, other):
        raise NotImplementedError

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "{}({})".format(self.__class__, self.name)

    def __str__(self):
        return str(self.name)

    def __invert__(self): 
        return self.complement

    def __add__(self): 
        raise NotImplementedError

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

class DL_Domain(ABC_Domain):
    """ An *instance* of a domain-level domain, and a *type* of sequence-level domain.

    This object supports two domain-level abstraction layers:
        - Layer 1: *short* and *long* (assuming "typical" sequences)
        - Layer 2: length (assuming "typical" sequences)
    ... for other layers see SL_Domain.

    Layer 1:
        A domain-level domain (DL_Domain) has two types: *short* and *long*.
        This destinction is motivated by domain-level DSD system design, where
        short domains (toehold domains) bind reversibly and long domains
        (branch migration domains) are bound irreversibly.  Based on these
        assumtions, a "rate-independent" domain-level reaction network can be
        enumerated. 

    Layer 2:
        A domain-level domain gets a *length*, which can be used to calculate 
        reaction rates assuming "typical" sequences, and, potentially, to 
        adjust the lengths for individual DL_Domain instances.

    Args:
        name(str): The name of a domain (must be unique!)
        dtype(str, optional): Specify whether domain is 'short' or 'long'. 
        length(int, optional): Specify domain length. 

    Raises:
        DSDDuplicationError: "Duplicate DL_Domain specification."
        DSDObjectsError: "Conflicting assignments of dtype and length."
        DSDObjectsError: "DL_Domain instance requires dtype and/or length"

    Class Variables:
        DTYPE_CUTOFF: assigns dtype from domain length
        SHORT_DOM_LEN: assigns length from dtype: 'short'
        LONG_DOM_LEN: assigns length from dtype: 'long'
        MEMORY: Used to initialize complementary domains and to raise Errors
            if domains are duplicated.

    """

    DTYPE_CUTOFF = 8  # largest value for short
    SHORT_DOM_LEN = 5
    LONG_DOM_LEN = 15

    MEMORY = dict()

    def __init__(self, name, dtype=None, length=None):
        if dtype not in set(['short', 'long', None]):
            raise DSDObjectsError('Cannot interpret dtype:', dtype)
        if length: assert isinstance(length, int)

        # Basic memory alerts.
        if name in DL_Domain.MEMORY:
            other = DL_Domain.MEMORY[name]
            error = DSDDuplicationError('Duplicate DL_Domain specification:', name, other.name)
            error.existing = other
            raise error
        else :
            DL_Domain.MEMORY[name] = self

        self._name = name

        if dtype and length :
            if not (dtype == 'short') == (length <= DL_Domain.DTYPE_CUTOFF):
                del DL_Domain.MEMORY[name]
                raise DSDObjectsError('Conflicting assignments of dtype and length')
        elif dtype :
            length = DL_Domain.SHORT_DOM_LEN if dtype == 'short' else DL_Domain.LONG_DOM_LEN
        elif length :
            dtype = 'short' if length <= DL_Domain.DTYPE_CUTOFF else 'long'
        else :
            del DL_Domain.MEMORY[name]
            raise DSDObjectsError('DL_Domain instance {} requires dtype and/or length'.format(name))

        self._length = length
        self._dtype = dtype
        self._complement = None

    @property
    def name(self): 
        return self._name

    @property
    def length(self):
        return self._length

    @property
    def is_complement(self):
        return self._name[-1] == '*'

    @property
    def complement(self):
        # If we initialize the complement, we need to know the class.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                raise DSDObjectsError('Complement has not been initialized.')
        return self._complement

    @property
    def dtype(self):
        return self._dtype

    def __repr__(self):
        return 'DL_Domain({}, {}, {})'.format(self._name, self._dtype, self._length)

    def __str__(self):
        return str(self._name)

    def __eq__(self, other):
        return self is other

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self._name < other._name

    def __gt__(self, other):
        return self._name > other._name

    def __le__(self, other):
        return self._name <= other._name

    def __ge__(self, other):
        return self._name >= other._name

    def __hash__(self):
        return hash(self.name)

class SL_Domain(ABC_Domain):
    """ An *instance* of a sequence-level domain.

    This object supports one layer of domain-level abstractions:
        - Layer 3: domains with context-dependent sequence modifications.
    ... for other layers see DL_Domain.

    Two instances of the same domain, should compete for the same complementary
    domain. Hoever, the competition depends on local context (neighboring
    domains, etc.) and might require sequence variations between domains on
    different strands.  SL_Domains (sequence-level domains) are
    context-dependent, sequence-optimized instances of DL_Domains (domain-level
    domains).  Complementarity and sequence-contraints are defined by the
    DL_Domain type.

    Class Variables:
        MEMORY: Used to initialize complementary domains and to raise Errors
            if domains are duplicated.
    """
    MEMORY = dict(dict())

    def __init__(self, dtype, sequence, variant=''):
        assert isinstance(dtype, DL_Domain)
        self._dtype = dtype
        self._sequence = self.enforce_constraints(sequence)

        if variant :
            self._name = dtype.name + '_' + variant
        else :
            self._name = dtype.name

        # Initialize MEMORY and alert duplicates.  Always do this at the end,
        # otherwise you might fill the memory with crap.
        if self._dtype.name in SL_Domain.MEMORY:
            if self._name in SL_Domain.MEMORY[self._dtype.name]:
                other = SL_Domain.MEMORY[self._dtype.name][self._name]
                error = DSDDuplicationError('Duplicate DL_Domain specification: {} {}'.format(self._name, other.name))
                error.existing = other
                raise error
            else :
                SL_Domain.MEMORY[self._dtype.name][self._name] = self
        else :
            SL_Domain.MEMORY[self._dtype.name] = {self._name : self}


    @property
    def name(self):
        return self._name

    @property
    def dtype(self):
        return self._dtype

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        raise NotImplementedError

    def enforce_constraints(self, constraint):
        """As for now, constraints are only on sequence length. """
        # make it strict for now:
        if not all(isinstance(c, str) for c in constraint):
            raise NotImplementedError('composite domains')
        if len(constraint) != len(self._dtype):
            raise DSDObjectsError('Sequence constraint violates domain type.', 
                    constraint, self._dtype)
        return constraint
    
    @property
    def complement(self):
        dtcomp = self._dtype.complement
        if dtcomp.name in SL_Domain.MEMORY:
            return list(SL_Domain.MEMORY[dtcomp.name].values())
        else :
            raise DSDObjectsError('Complement has not been initialized.')

    @property
    def variants(self):
        return list(SL_Domain.MEMORY[self._dtype.name].values())

    @property
    def length(self):
        return len(self._sequence)

    def is_logic(self, other):
        """ Domains are logically equivalent if their dtypes are equal.

        Note: This means that two instances of domains may be equal but differ
        on the sequence level. This is intended, as sequence-level modifications
        may be used to fine-tune context-dependent reaction dynamics.
        """
        return self._dtype is other._dtype

    def __repr__(self):
        return 'SL_Domain({}, {}, {})'.format(self._name, self._dtype, self._sequence)

    def __str__(self):
        return str(self._name)

    def __eq__(self, other):
        """ Domains are equal if their names are equal

        Note: This means that two instances of domains may be equal but differ
        on the sequence level. This is intended, as sequence-level modifications
        may be used to fine-tune context-dependent reaction dynamics.
        """
        if not isinstance(self, SL_Domain) or not isinstance(other, SL_Domain):
            return False
        return self.name == other.name

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self._name < other._name

    def __gt__(self, other):
        return self._name > other._name

    def __le__(self, other):
        return self._name <= other._name

    def __ge__(self, other):
        return self._name >= other._name

    def __hash__(self):
        return hash(self.name)

class DSD_Complex(object):
    """A sequence and structure pair.

    Other than domains, Complexes don't need the memory management, but it is
    useful in any case.  For example, if the same complex is produced twice in
    a different rotation, then this raises a DSDDuplicationError that contains
    both the exisiting complex object and the differences in rotations.

    Sequence and structure can be specified on the domain or on the nucleotide
    level, but they have to be of the same length. Although not implemented,
    one may define special characters for secondary structures that are more
    diverse than a regular dot-bracket string, e.g. 'h' = '(', '.', ')'

    Class Variables:
        ID: a counter to assing the next automatic name in a system.
        NAMES:
        MEMORY:

    Args:
        sequence (list): A domain-level or nucleotide-level sequence.
            TODO: enable to pass a DSD_StrandOrder instead.
        structure (list): A domain-level or nucleotide-level dot-bracket notation.
        name (str, optional): Name of this domain. If not specified, an automatic
            name is generated.
        prefix (str, optional): A prefix for automatic naming of Domains.
            Defaults to 'cplx'.
        memorycheck (bool, optional): Use built-in memory checks. Defaults to True.

    """

    ID = 0          # ID is used to assign names automatically
    NAMES = dict()  # NAMES[name] = canonical_form
    MEMORY = dict() # MEMORY[canonical_form] = self

    def __init__(self, sequence, structure, name='', prefix='cplx', memorycheck=True):
        # Assign name
        if name:
            self._name = name
        else:
            if prefix == '':
                raise DSDObjectsError('DSD_Complex prefix must not be empty!')
            if prefix[-1].isdigit():
                raise DSDObjectsError('DSD_Complex prefix must not end with a digit!')
            self._name = prefix + str(DSD_Complex.ID)
            DSD_Complex.ID += 1

        if len(sequence) != len(structure):
            raise DSDObjectsError(
                "DSD_Complex() sequence and structure must have same length")

        # should remain in same order as initialized
        self._sequence = sequence
        self._structure = structure
        self._canonical_form = None
        self._rotations = None
        self._strand_lengths = None

        # Initialized on demand:
        self._pair_table = None
        self._loop_index = None
        self._exterior_loops = None
        self._lol_sequence = None
        self._strands = None
        self._domains = None
        self._exterior_domains = None
        self._enclosed_domains = None

        # rotate strands into a canonical form.
        # two complexes are equal if they have equal canonial form
        self._memorycheck = memorycheck
        if self._memorycheck:
            canon = self.canonical_form # raises duplication error
            if self._name not in DSD_Complex.NAMES:
                DSD_Complex.NAMES[self._name] = canon
            else :
                raise DSDObjectsError('Duplicate DSD_Complex name!', self._name)
            DSD_Complex.MEMORY[self.canonical_form] = self

    def do_memorycheck(self, current = None, rotations=None):
        if current is None: 
            current = self.canonical_form 

        if current in DSD_Complex.MEMORY:
            other = DSD_Complex.MEMORY[current]
            error = DSDDuplicationError('Duplicate Complex specification:', 
                    current)
            error.existing = other
            error.rotations = abs(rotations - self.size) - other._rotations
            raise error

    @property
    def name(self):
        """str: name of the complex object. """
        return self._name

    @name.setter
    def name(self, name):
        if name in DSD_Complex.NAMES:
            raise DSDObjectsError('Duplicate DSD_Complex name!', name)
        if self._name in DSD_Complex.NAMES:
            del DSD_Complex.NAMES[self._name]
        self._name = name
        DSD_Complex.NAMES[self._name] = self.canonical_form

    @property
    def sequence(self):
        """list: sequence the complex object. """
        return self._sequence[:]

    @property
    def structure(self):
        """list: the complex structure. """
        return self._structure[:]

    @property
    def canonical_form(self):
        """ Sort sequence & structure in a unique way for each complex.

        Complexes are equal if they have the same canonical form:
            1) sort by sequence string
            2) sort by structure string

        Sets a variable self._rotations.
        """
        if not self._canonical_form:
            all_variants = dict()
            for e, new in enumerate(self.rotate(), 1):
                canon = tuple((tuple(map(str,self._sequence)), tuple(self._structure)))
                if canon not in all_variants:
                    all_variants[canon] = e
                    if self._memorycheck:
                        self.do_memorycheck(canon, e)
            self._canonical_form = sorted(list(all_variants.keys()), key = lambda x:(x[0],x[1]))[0]
            # Invert rotations to be compatible with 'wrap'
            self._rotations = abs(all_variants[self._canonical_form] - self.size)
        return self._canonical_form

    @property
    def size(self):
        """Size of the complex is the number of strands."""
        if not self._strand_lengths :
            if not self._lol_sequence:
                self._lol_sequence = utils.make_lol_sequence(self._sequence)
            self._strand_lengths = list(map(len, self._lol_sequence))
        return len(self._strand_lengths)
 
    @property
    def rotations(self):
        """Cyclic permutations between representation and canonical form."""
        return self._rotations

           
    def strand_length(self, pos):
        if not self._strand_lengths :
            if not self._lol_sequence:
                self._lol_sequence = utils.make_lol_sequence(self._sequence)
            self._strand_lengths = list(map(len, self._lol_sequence))
        return self._strand_lengths[pos]
 
    def rotate(self):
        """Generator function yields every rotation of the complex. """
        for i in range(self.size):
            yield self.rotate_once()

    def rotate_once(self):
        """Rotate the strands within the complex and return the updated object. """
        if "+" in self._sequence:
            p = self._sequence.index('+')
            self._sequence = self._sequence[p + 1:] + ["+"] + self._sequence[:p]

            tmpstruct = self.structure

            stack = []
            for i in range(p):
                if tmpstruct[i] == "(":
                    stack.append(i)
                elif tmpstruct[i] == ")":
                    try :
                        stack.pop()
                    except IndexError:
                        raise DSDObjectsError(
                                "Unbalanced parenthesis in secondary structure.")
            for i in stack:
                tmpstruct[i] = ")"

            stack = []
            for i in reversed(range(p + 1, len(tmpstruct))):
                if tmpstruct[i] == ")":
                    stack.append(i)
                elif tmpstruct[i] == "(":
                    try :
                        stack.pop()
                    except IndexError:
                        raise DSDObjectsError(
                                "Unbalanced parenthesis in secondary structure.")
            for i in stack:
                tmpstruct[i] = "("
            self._structure = tmpstruct[p + 1:] + ["+"] + tmpstruct[:p]
        self._pair_table = None
        self._loop_index = None
        self._lol_sequence = None
        self._exterior_domains = None
        return self

    def rotate_pairtable_loc(self, loc, n=None):
        """Maps the locus of a given pair-table to a new rotation. 

        For example, peppercorn finds a possible interaction between two
        complexes, but didn't guess the "correct" rotation. After DSD_Complex
        complains that the canonical form exists already, one can access the
        existing complex and the rotational distance from the error object. 

        """
        def wrap(x, m):
            """
            Mathematical modulo; wraps x so that 0 <= wrap(x,m) < m. x can be negative.
            """
            return (x % m + m) % m

        if n is None:
            n = self.size
        return (wrap(loc[0] + n, self.size), loc[1])

    @property
    def kernel_string(self):
        """str: print sequence and structure in `kernel` notation. """
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

    @property
    def domains(self):
        """ Returns a list of domains (sorted/unique) """
        if not self._domains:
            try :
                self._domains = set(self._sequence)
                if '+' in self._domains: self._domains.remove('+')
                self._domains = sorted(self._domains, key = lambda x: x.name)
            except AttributeError:
                raise DSDObjectsError('Not a domain-level sequence.', self.sequence)
        return self._domains

    @property
    def strands(self):
        raise NotImplementedError

    # Move set
    @property
    def pair_table(self):
        """ returns a structure in multistranded pair-table format. """
        # Make a new pair_table every time, it might get modified.
        return utils.make_pair_table(self.structure)


    @property
    def loop_index(self):
        """ returns the loop index of a structure. """
        if not self._pair_table:
            self._pair_table = utils.make_pair_table(self.structure)
        return utils.make_loop_index(self._pair_table)

    @property
    def lol_sequence(self):
        """ Returns sequence as a list of lists, without the '+' separator.
        
         Example:
          ['d1', 'd2', '+', 'd3'] ==> [['d1', 'd2'], ['d3']]
        """
        return utils.make_lol_sequence(self._sequence)

    def get_loop_index(self, loc):
        if not self._pair_table:
            self._pair_table = utils.make_pair_table(self.structure)
        if not self._loop_index:
            self._loop_index, self._exterior_loops = utils.make_loop_index(self._pair_table)
        return self._loop_index[loc[0]][loc[1]]

    def get_domain(self, loc):
        if not self._lol_sequence:
            self._lol_sequence = utils.make_lol_sequence(self._sequence)
        return self._lol_sequence[loc[0]][loc[1]]
    
    def get_paired_loc(self, loc):
        """ 
        Returns the paired element in the pair-table. 
        Raises: IndexError if there are negative elements in loc
        """
        if loc[0] < 0 or loc[1] < 0:
            raise IndexError
        if not self._pair_table:
            self._pair_table = utils.make_pair_table(self.structure)
        return self._pair_table[loc[0]][loc[1]]

    @property
    def enclosed_domains(self):
        if not self._enclosed_domains:
            # Ask for exterior domains, that will also
            # fill the enclosed domains
            _ = self.exterior_domains
        return self._enclosed_domains

    @property
    def exterior_domains(self):
        """
        Returns all domains in exterior loops.
        """
        if not self._exterior_domains:
            if not self._pair_table:
                self._pair_table = utils.make_pair_table(self.structure)
            if not self._loop_index or self._exterior_loops:
                self._loop_index, self._exterior_loops = utils.make_loop_index(self._pair_table)

            self._exterior_domains = []
            self._enclosed_domains = []
            for si, strand in enumerate(self._loop_index):
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

        TODO: breaks if used on strings, due to the inverse operator
        """
        if not self._pair_table:
            self._pair_table = utils.make_pair_table(self.structure)
        for si, strand in enumerate(self._pair_table):
            for di, domain in enumerate(strand):
                loc = (si,di)
                cloc = self._pair_table[si][di] 
                if not (cloc is None or self.get_domain(loc) == ~self.get_domain(cloc)):
                    raise DSDObjectsError('Domains cannot pair', 
                            self.get_domain(loc), self.get_domain(cloc))
        return True

    @property
    def is_connected(self):
        if not self._pair_table:
            self._pair_table = utils.make_pair_table(self.structure)
        if not self._loop_index:
            try :
                # make loop index raises error when disconnected.
                self._loop_index, self._exterior_loops = utils.make_loop_index(self._pair_table)
            except utils.DSDUtilityError as e:
                return False
        return True

    def __repr__(self):
        return 'DSD_Complex({}, {})'.format(self.name, self.canonical_form)

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        """ Test if two complexes are equal.

        They are equal if they have the same sequence and secondary structure.
        The condition does not depend on MEMORY and can compare instances of
        different child-classes.
        
        """
        if not isinstance(self, DSD_Complex) or not isinstance(other, DSD_Complex):
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

    def __len__(self):
        # ambiguos... length of sequence? length of nucleotides?
        raise NotImplementedError

    def __hash__(self):
        return hash(self.canonical_form)

    def __add__(self): 
        raise NotImplementedError

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)


class DSD_Macrostate(object):
    """
    A set of complexes.
    """

    NAMES = dict()  # NAMES[name] = complex.canonical_form
    MEMORY = dict() # MEMORY[complex.canonical_form] = self

    def __init__(self, complexes, name='', prefix='r', 
            representative=None, memorycheck=True):

        # Resting sets store complexes in canonical form?
        assert len(set(complexes)) == len(complexes)
        self._complexes = sorted(complexes, key=lambda x : x.canonical_form)
        
        if representative is None:
            self._canonical = self._complexes[0]
        else:
            assert isinstance(representative, DSD_Complex)
            assert representative in set(complexes)
            self._canonical = representative

        # Assign name
        if name:
            self._name = name
        else:
            self._name = prefix + self._canonical.name

        # two complexes are equal if they have equal canonial form
        if memorycheck :
            if self.canonical_form in DSD_Macrostate.MEMORY:
                other = DSD_Macrostate.MEMORY[self.canonical_form]
                if other != self :
                    error = DSDObjectsError(
                            'Conflicting Macrostate Assignment: {} and {}.'.format(
                                self, other))
                    error.existing = other
                    raise error
                else :
                    error = DSDDuplicationError('Duplicate Macrostate specification:', 
                            self, other)
                    error.existing = other
                    raise error
            elif self._name in DSD_Macrostate.NAMES:
                raise DSDObjectsError('Duplicate Macrostate name!', self._name)
            else :
                for canon in map(lambda x: x.canonical_form, self._complexes):
                    DSD_Macrostate.MEMORY[canon] = self
                DSD_Macrostate.NAMES[self._name] = self.canonical_form


    @property
    def name(self):
        return self._name

    @property
    def complexes(self):
        """
        A list of complexes in the resting set
        """
        return self._complexes[:]

    @property
    def canonical_name(self):
        """
        Gives the canonical name of the resting set, chosen by its "first" 
        element. (sorted by canonical forms)
        """
        return self._canonical.name

    @property
    def canonical(self):
        """
        See ``canonical_name``.
        """
        return self._canonical

    @property
    def canonical_form(self):
        """
        See ``canonical_name``.
        """
        return self._canonical.canonical_form

    @property
    def kernel_string(self):
        return self.canonical.kernel_string

    def __repr__(self):
        return "DSD_Macrostate(\"%s\", %s)" % (self.name, str(self.complexes))

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        """
        Two resting sets are equal if their complexes are equal
        """
        return (self.complexes == other.complexes)

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

    def __len__(self):
        """ The number of species in a macrostate """
        return len(self._complexes)

    def __hash__(self):
        return hash(self.canonical_form)

    def __add__(self): 
        raise NotImplementedError

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

class DSD_Reaction(object):
    """ A reaction pathway.

    Args:
      reactants (list): A list of reactants. Reactants can be 
        :obj:`DSD_Macrostate()` or :obj:`DSD_Complex()` objects.
      products (list): A list of products. Products can be
        :obj:`DSD_Macrostate()` or :obj:`DSD_Complex()` objects.
      rtype (str, optional): Reaction type, e.g. bind21, condensed, .. Defaults to None.
      rate (flt, optional): Reaction rate. A reaction rate 

    TODO: think about reversible reactions (invert operator, etc)
    """

    MEMORY = dict()
    Rate = namedtuple('rate', 'constant units')

    def __init__(self, reactants, products, rtype=None, rate=None, memorycheck=True):
        if not(all(isinstance(c, 
            (DSD_Complex, DSD_Macrostate)) for c in reactants + products)):
            raise DSDObjectsError('Reactants and products must be DSD Objects')
        self._reactants = reactants
        self._products = products
        self._rtype = rtype

        if rate is None or isinstance(rate, DSD_Reaction.Rate):
            self._rate = rate 
        else: 
            self._rate = DSD_Reaction.Rate(rate, ['M'] * (self.arity[0] - 1) + ['s']) 
        
        # Used for __eq__ 
        self._canonical_form = None

        if memorycheck:
            if self.canonical_form in DSD_Reaction.MEMORY:
                other = DSD_Reaction.MEMORY[self.canonical_form]
                error = DSDDuplicationError('Duplicate reaction in system.', self.canonical_form)
                error.existing = other
                raise error
            else :
                DSD_Reaction.MEMORY[self.canonical_form] = self

    @property
    def rate(self):
        """Return reaction rate constant and units."""
        return self._rate

    @rate.setter
    def rate(self, value):
        """Set reaction rate constant and units."""
        assert isinstance(value, DSD_Reaction.Rate)
        self._rate = value

    @property
    def const(self):
        """Return reaction rate constant."""
        return self._rate if self._rate is None else self._rate.constant

    @const.setter
    def const(self, value):
        """Set reaction rate with SI units."""
        self._rate = DSD_Reaction.Rate(value, ['M'] * (self.arity[0] - 1) + ['s']) 

    @property
    def rateunits(self):
        """Return reaction rate units."""
        return ''.join(map('/{}'.format, self._rate.units))

    def rateformat(self, output_units):
        """Set reaction rate constant and units."""
        assert isinstance(output_units, list)
        newc = self.const
        for i,o in zip(self._rate.units, output_units):
            newc = utils.convert_units(newc, o, i) # 1/M 1/s
        return DSD_Reaction.Rate(newc, output_units)

    @property
    def rtype(self):
        """str: *peppercorn* reaction type (bind21, condensed, ...) """
        return self._rtype

    @property
    def reactants(self):
        """list: list of reactants. """
        return self._reactants[:]

    @property
    def products(self):
        """list: list of products. """
        return self._products[:]

    @property
    def arity(self):
        """(int, int): number of reactants, number of products."""
        return (len(self._reactants), len(self._products))

    @property
    def kernel_string(self):
        return "{} -> {}".format("  +  ".join(r.kernel_string for r in self.reactants),
                "  +  ".join(p.kernel_string for p in self.products))

    @property
    def canonical_form(self):
        if not self._canonical_form:
            react = tuple(sorted(list(map(lambda x: x.canonical_form, self.reactants))))
            prods = tuple(sorted(list(map(lambda x: x.canonical_form, self.products))))

            self._canonical_form = tuple((react, prods, self.rtype))
        return self._canonical_form

    def __repr__(self):
        return "{}({})".format(self.__class__, self.canonical_form)

    def __str__(self):
        """prints the formal chemical reaction."""
        return "{} -> {}".format(
            " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

    def __eq__(self, other):
        """bool: Checks if DSD_Reaction objects have the same rtype, reactants,
        and products. rate can be different.""" 
        if not isinstance(self, DSD_Reaction) or not isinstance(other, DSD_Reaction):
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

class DSD_StrandOrder(object):
    """ A single strand, or a list of strands in specific order.

    Cyclic permutations of a strand order are equivalent.

    The sequence of a strand order must be a sequence of domains.

    """

    ID = 0          # ID is used to assign names automatically
    NAMES = dict()  # NAMES[name] = canonical_form
    MEMORY = dict() # MEMORY[canonical_form] = self

    # Canonical form for strand order
    def __init__(self, sequence, name='', prefix='StrandOrder_', memorycheck=True):

        # Assign name
        if name:
            self._name = name
        else:
            if prefix == '':
                raise DSDObjectsError('DSD_StrandOrder prefix must not be empty!')
            if prefix[-1].isdigit():
                raise DSDObjectsError('DSD_StrandOrder prefix must not end with a digit!')
            self._name = prefix + str(DSD_StrandOrder.ID)
            DSD_StrandOrder.ID += 1

        # should remain in same order as initialized
        self._sequence = sequence
        self._canonical_form = None
        self._rotations = None
        self._lol_sequence = None

        # Initialized on demand:
        self._strand_lengths = None
        self._strands = None
        self._domains = None

        # rotate strands into a canonical form.
        # two complexes are equal if they have equal canonial form
        self._memorycheck = memorycheck
        if self._memorycheck:
            canon = self.canonical_form # raises duplication error
            if self._name not in DSD_StrandOrder.NAMES:
                DSD_StrandOrder.NAMES[self._name] = canon
            else :
                raise DSDObjectsError('Duplicate DSD_StrandOrder name!', self._name)
            DSD_StrandOrder.MEMORY[self.canonical_form] = self

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if name in DSD_StrandOrder.NAMES:
            raise DSDObjectsError('Duplicate DSD_StrandOrder name!', name)
        if self._name in DSD_StrandOrder.NAMES:
            del DSD_StrandOrder.NAMES[self._name]
        self._name = name
        DSD_StrandOrder.NAMES[self._name] = self.canonical_form

    @property
    def sequence(self):
        """list: sequence the complex object. """
        return self._sequence[:]

    @property
    def kernel_string(self):
        """str: print sequence and structure in `kernel` notation. """
        return "{}".format(" ".join(map(str,self._sequence)))

    @property
    def size(self):
        """Size of the StrandOrder is the number of strands."""
        return len(self._lol_sequence)

    @property
    def canonical_form(self):
        """ Sort Sequences in a unique way for each complex.
        """
        if not self._canonical_form:
            all_variants = dict()
            if not self._lol_sequence:
                self._lol_sequence = utils.make_lol_sequence(self._sequence)
            for e in range(0, self.size):
                canon = tuple(map(lambda x: tuple(map(str, x)), 
                    self._lol_sequence[e:] + self._lol_sequence[:e]))
                if canon not in all_variants:
                    all_variants[canon] = e
                    if self._memorycheck:
                        self.do_memorycheck(canon, e)
            self._canonical_form = sorted(all_variants)[0]
            self._rotations = abs(all_variants[self._canonical_form] - self.size)
        return self._canonical_form
 
    @property
    def rotations(self):
        """Cyclic permutations between representation and canonical form."""
        return self._rotations

    @property
    def lol_sequence(self):
        """ Returns sequence as a list of lists, without the '+' separator.
        
         Example:
          ['d1', 'd2', '+', 'd3'] ==> [['d1', 'd2'], ['d3']]
        """
        return utils.make_lol_sequence(self._sequence)

    def do_memorycheck(self, current = None, rotations=None):
        if current is None: 
            current = self.canonical_form 

        if current in DSD_StrandOrder.MEMORY:
            other = DSD_StrandOrder.MEMORY[current]
            error = DSDDuplicationError('Duplicate DSD_StrandOrder specification:', 
                    current)
            error.existing = other
            error.rotations = abs(rotations - self.size) - other._rotations
            raise error

    def __repr__(self):
        return 'DSD_StrandOrder({}: {})'.format(self.name, ' '.join(map(str, self.sequence)))

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
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

    def __len__(self):
        # ambiguos... length of sequence? number of strands?
        raise NotImplementedError

    def __hash__(self):
        return hash(self.canonical_form)

# Deprecated
class DSD_RestingSet(DSD_Macrostate):
    def __init__(self, *kargs, **kwargs):
        print("DEPRECATION WARNING: DSD_RestingSet is now called DSD_Macrostate.")
        super(DSD_RestingSet, self).__init__(*kargs, **kwargs)

