#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name='dsdobjects',
    version='0.7.1',
    description='Base classes for DSD design',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url='https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects',
    author='Stefan Badelt',
    author_email='badelt@caltech.edu',
    license='MIT',
    download_url = 'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects/archive/v0.7.1.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 2.7',
        ],
    install_requires=['future'],
    packages=['dsdobjects', 'dsdobjects.parser', 'dsdobjects.core'],
    test_suite='tests',
)

