#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

# Dynamically figure out the version
setup(
    name='dsdobjects',
    version='0.6',
    description='Base classes for DSD design',
    long_description=readme,
    url='https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects',
    author='Stefan Badelt',
    author_email='badelt@caltech.edu',
    license=license,
    download_url = 'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects/archive/v0.6.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        ],
    packages=['dsdobjects', 'dsdobjects.parser'],
    test_suite='tests',
)

