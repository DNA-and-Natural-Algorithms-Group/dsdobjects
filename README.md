# dsdobjects: an object library for DSD programming

This Python module provides a library of protype objects and base classes for
domain-level strand displacement (DSD) programming. There are two types of
usage: 
 1) ready-to-go prototype objects, 
 2) tweak-em-yourself core objects.

If you start with a new project, you want to use the prototypes. These are
out-of-the box functional classes that can be initialized, for example, using
various standards of the text file format *.pil. 

## Installation
To install this library use pip:
```
$ pip install dsdobjects
```
or the following command in the root directory:
```
$ python ./setup.py install
```


### Quick Start with object prototypes.
```py
from dsdobjects import SequenceConstraint, StrandOrder, LogicDomain, Domain, Complex, Macrostate, Reaction
```

```py
 # Define a few toy domains:
 a = LogicDomain('a', dtype='long')
 b = LogicDomain('b', dtype='long', length=9)
 t = LogicDomain('t', dtype='short', length=6)
 
 # LogicDomains have exactly one complement, it can be initialized 
 # and/or accessed using the __invert__ operator. The built-in 
 # memory management ensures that there is only one object for each domain.
 assert (a is ~(~a))
 
 # Use the Domains to define a Complex ...
 cplx = Complex([a, b, c, ~b, '+', ~a], list('((.)+)'), name='rudolf')

 # ... and test some of the built-in complex properties:
 cplx.kernel_string
 cplx.canonical_form
 cplx.size
 for r in cplx.rotate():
     print r.kernel_string
 cplx.pair_table
 
 # If you were to define two complexes as one disconnected complex ... 
 cplx = Complex([a, b, c, ~b, '+', ~a], list('.(.)+.'), name='cplx')
 assert cplx.is_connected is False

 # ... you can quickly and return the indiviudal complexes:
 cx1, cx2 = cplx.split()
```

### Quick Start with text input.
For example load an output file from peppercorn into dsdobjects:

```py
 from dsdobjects import read_pil
 
 domains, complexes, macrostates, detailed_rxns, condensed_rxns = read_pil(filename.pil)
 
 myobject = read_pil_line("length d5 = 7")
 assert isinstance(myobject, LogicDomain)
 
 myobject = read_pil_line("sequence d5 = NNNNN")
 assert isinstance(myobject, Domain)
```

## Abut the core objects.
If prototypes are not sufficient, you can make your own objects by inheriting
from the core objects. Core objects provide a basic set of __builtin__
functions (e.g. equality, sorting), basic properties, memory management.  One
way to get started is by copying the prototypes file into your project, and
adapt it to your needs. Consider a pull request back into the main dsdobjects
repository!


### Quick Start with core objects.
Inheritance from dsdobjects.base_classes provides only basic functions such as
'~', '==', '!=', and access to the built-in memory management for each class.
Some potential ambiguities, such as requesting the complement of a Domain,  or
the length of a complex must be defined upon inheritance.

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

## Version

0.7 -- Python 3.x support / prototypes
  * basic support of prototype objects
  * added StrandOrder base_class and prototpye
  * allow parsing of infinite error bars for reaction rates
  * DSD_Restingset renamed to DSD_Macrostate
  * broken backward compatibility:
      reaction rates are now namedtuples

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

