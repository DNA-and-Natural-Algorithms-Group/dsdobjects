#
# dsdobjects/prototypes.py
#   - copy and/or modify together with tests/test_prototypes.py
#
import logging
log = logging.getLogger(__name__)

from .singleton import SingletonError
from .base_classes import ObjectInitError, DomainS, ComplexS
from .complex_utils import rotate_complex_db, make_strand_table, wrap

class DomainA(DomainS):
    # Allows for automatic name assignments.
    DTYPE_CUTOFF = 8 
    SHORT_DOM_LEN = 5
    LONG_DOM_LEN = 15
    ID = 1

    @classmethod
    def identifiers(cls, name = None, length = None, prefix = 'd', dtype = None):
        my = cls
        """ tuple: A method that must be accessible without initializing the object. """
        if name is None:
            name = f'{prefix}{my.ID}'
        if length is None:
            length = my.SHORT_DOM_LEN if dtype == 'short' else my.LONG_DOM_LEN if dtype == 'long' else None
        elif dtype and not (dtype == 'short') == (length <= my.DTYPE_CUTOFF):
            raise ObjectInitError(f'Conflicting arguments {dtype} and {length}.')
        if length is None and name[-1] == '*':
            try:
                length = len(cls(name[:-1]))
            except SingletonError:
                pass
        return ((name, length), name, {}) if length is not None else (None, name, {})

    def __init__(self, name = None, length = None, prefix = 'd', dtype = None):
        my = self.__class__ # same as in identifiers!
        if name is None:
            name = f'{prefix}{my.ID}'
            my.ID += 1
        if length is None:
            length = my.SHORT_DOM_LEN if dtype == 'short' else my.LONG_DOM_LEN if dtype == 'long' else None
        super(DomainA, self).__init__(name, length)

    @property
    def dtype(self):
        return 'short' if self.length <= self.__class__.DTYPE_CUTOFF else 'long'

    @property
    def complement(self):
        """ obj: the complementary domain object. """
        return self.__class__(self.cname, self.length)

class ComplexA(ComplexS):
    # Allows for automatic name assignments.
    ID = 1

    @classmethod
    def identifiers(cls, sequence, structure, name = None, prefix = 'cplx', **kwargs):
        """ tuple: A method that must be accessible without initializing the object. """
        my = cls
        if sequence is None or structure is None:
            if name is None:
                raise ObjectInitError('Insufficient arguments for Complex initialization.')
            return (None, name, {})
        if name is None:
            name = f'{prefix}{my.ID}'
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

    def __init__(self, sequence, structure, name = None, prefix = 'cplx', **kwargs):
        my = self.__class__
        if name is None:
            name = f'{prefix}{my.ID}'
            my.ID += 1
        super(ComplexA, self).__init__(sequence, structure, name = name, **kwargs)

