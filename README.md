# dsdobjects - classes / prototypes for DSD design
This module provides standardized Python classes for domain-level strand
displacement programming:

- SequenceConstraint
- DL_Domain
- SL_Domain
- DSD_Complex
- DSD_Reaction
- DSD_Macrostate
- DSD_StrandOrder

Using the available prototype classes provides a quick but standardized way to
write new algorithms for DSD programming. Alternatively, inheritance from
dsdobjects.base_classes provides only basic functions such as '~', '==', '!=',
and access to the built-in memory management for each class. Some potential
ambiguities, such as requesting the complement of a Domain,  or the length of a
complex must be defined upon inheritance.

## Quick Start with object prototypes.

```py
from dsdobjects import DL_Domain as LogicDomain
from dsdobjects import DSD_Complex as Complex

# Define a few toy domains:
a = LogicDomain('a', dtype='long')
b = LogicDomain('b', dtype='long', length=9)
t = LogicDomain('t', dtype='short', length=6)

# DL_Domains have always only one complement, it can be 
# initialized and/or accessed using the __invert__ operator.
assert (a is ~(~a))


# Use the Domains to define a Complex:
>>> cplx = Complex([a, b, c, ~b, '+', ~a], list('((.)+)'), name='rudolf')
cplx.kernel_string
cplx.canonical_form
cplx.size
for r in cplx.rotate():
    print r.kernel_string
cplx.pair_table

# Define a two disconnected Complexes... 
cplx = Complex([a, b, c, ~b, '+', ~a], list('.(.)+.'), name='cplx')
cx1, cx2 = cplx.split()


from dsdobjects import SequenceConstraint as Sequence
from dsdobjects import SL_Domain as SequenceDomain

seq1 = Sequence('ACGTNNGT', molecule='DNA')
seq2 = Sequence('HHHHHHHH', molecule='DNA')
seq3 = seq1 + seq2
seq1.add_constraint('HHHHHHHH')
print seq1.constraint
print seq1.complement
print seq1.wc_complement
print seq1.reverse_complement
print seq1.reverse_wc_complement

seq1c = seq1.wc_complement
```


```py
from dsdobjects import DL_Domain

# A personalized domain that extends the DL_Domain base class.
class MyDomain(DL_Domain):

    def __init__(self, name, dtype=None, length=None):
        super(MyDomain, self).__init__(name, dtype, length)
 
    @property
    def complement(self):
        # Automatically initialize or return the complementary domain.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = MyDomain(cname, self.dtype, self.length)
        return self._complement

```

Inheriting from the DL_Domain base class enables standardized built in
functions such as '~', '==', '!=', and provides a built-in memory management
raising the DSDDuplicationError when conflicting domain names are chosen.


```py
>>> # Initialize a Domain.
>>> x = MyDomain('hello', dtype='short')
>>> # The '~' operator calls x.complement
>>> y = ~x
>>> (y == ~x)
True

```

These and many more functionalities and sanity checks are also available for
other objects. See the respective docstrings for details.  

## Install
To install this library, use the following command in the root directory:
```
$ python ./setup.py install
```
or use local installation:
```
$ python ./setup.py install --user
```

## Version
0.7 -- Python 3.x support / prototypes
    * basic support of prototype objects
    * added StrandOrder base_class and prototpye
    * allow parsing of infinite error bars for reaction rates
    * DSD_Restingset renamed to DSD_Macrostate
    * broken backward compatibility:
        reaction rates are now named tuples

0.6.3 -- added parser for seesaw language

0.6.2 -- bugfix for restingsets with given representative

0.6.1 -- adapted setup.py when used as pypi dependency

0.6 -- PIL parser supports concentration format
  * "non-equal" bugfixes in base_classes.py
  * supports rate-error bars when parsing PIL format

0.5 -- improved canonical forms

## Author
Stefan Badelt

### Contributors
This library contains adapted code from various related Python packages coded
in the [DNA and Natural Algorithms Group], Caltech:
  * "DNAObjects" coded by Joseph Berleant and Joseph Schaeffer 
  * [peppercornenumerator] coded by Kathrik Sarma, Casey Grun and Erik Winfree
  * [nuskell] coded by Seung Woo Shin

## Projects depending on dsdobjects
  * [peppercornenumerator]
  * [nuskell]


## License
MIT

[nuskell]: <http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell>
[peppercornenumerator]: <http://www.github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[DNA and Natural Algorithms Group]: <http://dna.caltech.edu>
