# dsdobjects: an object library for DSD programming

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/dna-and-natural-algorithms-group/dsdobjects)](https://github.com/dna-and-natural-algorithms-group/dsdobjects/tags)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/dna-and-natural-algorithms-group/dsdobjects?include_prereleases)](https://github.com/dna-and-natural-algorithms-group/dsdobjects/releases)
[![PyPI version](https://badge.fury.io/py/dsdobjects.svg)](https://badge.fury.io/py/dsdobjects)
[![PyPI - License](https://img.shields.io/pypi/l/dsdobjects)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/dna-and-natural-algorithms-group/dsdobjects.svg?branch=development)](https://travis-ci.com/github/dna-and-natural-algorithms-group/dsdobjects)
[![Codecov branch](https://img.shields.io/codecov/c/github/dna-and-natural-algorithms-group/dsdobjects/development)](https://codecov.io/gh/dna-and-natural-algorithms-group/dsdobjects)

A library of base classes for domain-level strand displacement (DSD)
programming.  If you are starting a new project that requires domains,
complexes, secondary structures, nucleotide sequences, reactions, etc.  then
you might find this module useful. Note that this is the in-house library of
the [nuskell] compiler framework; it handels the parsing of supported input file
formats (mainly \*.PIL) and their translation into the respective DSD objects.

All objects provided here are **singletons**. In other words, once you
intialize a Domain ``d1 = DomainS(name = 'a', length = 15)``, you will not be
allowed to initialize a new domain with the same name but different length,
unless you delete (all references to) the variable ``d1`` first. If you try to
initialize the same domain a second time, then the new object ``is`` d1.
Generally, each (DomainS, ComplexS, MacrostateS, ReactionS) must be initialized
with parameters that define their *canonical form* and (optionally) a name. If
both are given, the library checks that there are no conflicts with existing
objects, if the name is not provided, the library tries to initialize a new
object with an automatic name, but raises a *SingletonError* which holds a
reference to the existing object. If only the name is given, the existing
object is returned.

This library is expected to evolve further, potentially breaking backward
compatibility as new challenges are waiting in the [nuskell] compiler
framework.  Don't hesitate to contact the authors with questions about future
plans. Inheritance from the provided objects is fully supported and encouraged.

## Installation
To install this library use pip:
```
$ pip install dsdobjects
```
or the following command in the root directory:
```
$ python ./setup.py install
```

### Quick Start
```py
from dsdobjects import DomainS, ComplexS

# Define a few toy domains:
a = DomainS('a', length = 15)
b = DomainS('b', length = 9)
c = DomainS('c', length = 6)

# DomainS objects have exactly one complement, it can be initialized 
# and/or accessed using the __invert__ operator. The singleton type
# ensures that there is only one object for each domain.
assert (a is ~(~a))

# Use the Domains to define a Complex ...
foo = ComplexS([a, b, c, ~b, '+', ~a], list('((.)+)'), name = 'foo')

# ... and test some of the built-in complex properties:
foo.kernel_string
foo.canonical_form
foo.size
foo.pair_table
for (se, ss) in foo.rotate():
    print(se, ss)

# If you initialize a disconnected complex ... 
bar = ComplexS([a, b, c, ~b, '+', ~a], list('.(.)+.'), name = 'bar')
assert bar.is_connected is False
# ... use split to get all indiviudal complexes:
cx1, cx2 = bar.split()
```

### Quick Start from PIL files
Initialize prototype objects by loading a system (or a single line) of \*.PIL
file format:

```py
from dsdobjects import DomainS
from dsdobjects.objectio import set_io_objects, clear_io_objects, read_pil, read_pil_line

# Use the builtin singleton obects form the dsdobjects library.
set_io_objects()

# The following dictionary contains references to all objects.
outdict = read_pil('filename', is_file = True)

# The following let's you quickly initialize single objects.
d5 = read_pil_line("length d5 = 7")
assert isinstance(d5, DomainS)

d6 = read_pil_line("sequence d6 = NNNNN")
assert isinstance(d6, DomainS)
assert d6.sequence == 'NNNNN'
```

## Version
0.8 -- requires Python<=3.7
  * complete rewrite of the library to use singleton objects with weakref

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

