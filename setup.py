#!/usr/bin/env python3

from setuptools import setup
import os
from subprocess import check_output

prefix = os.environ['HOME']

bindir = prefix + '/bin'
if not os.path.exists(bindir):
    os.makedirs(bindir)

check_output('make PREFIX="' + prefix + '" install', shell=True)

setup(
    name = 'sarks',
    version = '0.1',
    description = 'Suffix Array Kernel Smoothing',
    author = 'Dennis Wylie',
    author_email = 'denniswylie@austin.utexas.edu',
    py_modules = ['sarks'],
    install_requires = ['biopython', 'intervaltree', 'numpy', 'pandas', 'pyfaidx', 'scipy']
)
