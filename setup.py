#!/usr/bin/env python

from setuptools import setup, find_packages

LONG_DESCRIPTION="""
This module provides Python parent classes for domain-level strand displacement
programming, e.g:
 - SequenceConstraint
 - DL_Domain
 - SL_Domain
 - DSD_Complex
 - DSD_Reaction
 - DSD_RestingSet
"""

setup(
    name='dsdobjects',
    version='0.6.3',
    description='Base classes for DSD design',
    long_description=LONG_DESCRIPTION,
    url='https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects',
    author='Stefan Badelt',
    author_email='badelt@caltech.edu',
    license='MIT',
    download_url = 'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects/archive/v0.6.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        ],
    install_requires=['future'],
    packages=['dsdobjects', 'dsdobjects.parser'],
    test_suite='tests',
)

