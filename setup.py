#!/usr/bin/env python

# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

# Dynamically figure out the version
setup(
    name='dsdobjects',
    version='0.6',
    description='base classes for DSD design',
    long_description=readme,
    author='Stefan Badelt',
    author_email='badelt@caltech.edu',
    license=license,
    packages=['dsdobjects', 'dsdobjects.parser'],
    test_suite='tests',
)

