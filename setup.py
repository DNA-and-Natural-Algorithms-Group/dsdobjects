#!/usr/bin/env python

from setuptools import setup, find_packages

LONG_DESCRIPTION="""
This Python module provides a library of protype objects and base classes for
domain-level strand displacement (DSD) programming. There are two types of
usage: 1) ready-to-go prototype objects, 2) tweak-em-yourself core objects.

from dsdobjects import SequenceConstraint, LogicDomain, Domain, Complex, Reaction, Macrostate, StrandOrder
"""

setup(
    name='dsdobjects',
    version='0.7',
    description='Base classes for DSD design',
    long_description=LONG_DESCRIPTION,
    url='https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects',
    author='Stefan Badelt',
    author_email='badelt@caltech.edu',
    license='MIT',
    download_url = 'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects/archive/v0.7.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 2.7',
        ],
    install_requires=['future'],
    packages=['dsdobjects', 'dsdobjects.parser', 'dsdobjects.core'],
    test_suite='tests',
)

