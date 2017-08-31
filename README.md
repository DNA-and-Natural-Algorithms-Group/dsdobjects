# dsdobjects - base classes for DSD design

Programming domain-level strand displacement (DSD) systems can require
different notions of "domains", "complexes", "reaction networks", etc. This
module povides minimal standardized Python parent classes:

  - SequenceConstraint
  - DSD_DL_Domain
  - DSD_SL_Domain
  - DSD_Complex
  - DSD_Reaction
  - DSD_RestingSet
  - DSD_TestTube

These base classes are currently used in the projects *Nuskell* and
*Peppercornenumerator*. An example for extending a base class is shown below:

```py
from dsdobjects import DSD_DL_Domain

# A personalized domain that extends the DSD_DL_Domain base class, allowing for 
# easy initialization of complementary domains.
class MyDomain(DSD_DL_Domain):

    def __init__(self, name, dtype=None, length=None):
        super(MyDomain, self).__init__(name, dtype, length)
 
    @property
    def complement(self):
        # If we initialize the complement, we need to know the class.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = DummyDomain(cname, self.dtype, self.length)
        return self._complement

```

The above class MyDomain() uses standardized built in functions such as 

```py
>>> # Initialize a Domain.
>>> x = MyDomain('hello', dtype='short')
>>> # The '~' operator calls x.complement
>>> y = ~x
>>> (y == ~x)
True

```

All of these obects come with sanity checks, ensuring that only one instance
of a particular domain, complex or reaction is initialized. Sanity checks also 
ensure that sequence constraints of (complementary) DSD_SL_Domains are 
properly enforced, etc.

## Install

```
$ python ./setup.py install --user
```

## Author
Stefan Badelt

## License
MIT
