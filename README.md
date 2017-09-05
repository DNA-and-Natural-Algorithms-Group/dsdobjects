# dsdobjects - base classes for DSD design

This module provides minimal standardized Python parent classes for 
domain-level strand displacement programming:

  - SequenceConstraint
  - DL_Domain
  - SL_Domain
  - DSD_Complex
  - DSD_Reaction
  - DSD_RestingSet

These base classes are currently used in the projects [nuskell] and
[peppercornenumerator]. An example for extending a base class is shown below:

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
                self._complement = DummyDomain(cname, self.dtype, self.length)
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
0.1 -- initial release

## Author
Stefan Badelt

### Contributors
This library contains adapted code from various related Python packages coded
in the [DNA and Natural Algorithms Group], Caltech:
  * "DNAObjecs" coded by Joseph Berleant and Joseph Schaeffer 
  * [peppercornenumerator] coded by Kathrik Sarma, Casey Grun and Erik Winfree
  * [nuskell] coded by Seung Woo Shin


## License
MIT

[nuskell]: <http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell>
[peppercornenumerator]: <http://www.github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[DNA and Natural Algorithms Group]: <http://dna.caltech.edu>
