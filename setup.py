#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'dsdobjects',
    version = '0.8',
    description = 'Base classes and prototype objects for DSD design',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = "text/markdown",
    author = 'Stefan Badelt',
    author_email = 'bad-ants-fleet@posteo.eu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects',
    license = 'MIT',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.7',
    install_requires = ['pyparsing'],
    packages = find_packages(),
    test_suite = 'tests',
)

