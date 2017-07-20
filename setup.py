from setuptools import setup
import os
import site
from subprocess import check_output

## -----------------------------------------------------------------------
## prefix + '/bin' will be directory into which cpp executables are placed
## -----------------------------------------------------------------------
prefix = site.USER_BASE

bindir = prefix + '/bin'
if not os.path.exists(bindir):
    os.makedirs(bindir)

## -----------------------------------------------------------------
## incdir is passed to makefile as flag for suffix-array compilation
## to specify location of local seqan installation if necessary
## -----------------------------------------------------------------
# incdir = '-I/work/03319/dwylie/seqan/seqan-library-1.4.2/include'
incdir = ''

## -----------------------------------------------------------------------
## build cpp executables and move them to prefix + '/bin' by invoking make
## -----------------------------------------------------------------------
check_output('make INC="' + incdir + '" PREFIX="' + prefix + '" install', shell=True)

setup(
    name = 'sarks',
    version = '0.1',
    description = 'Suffix Array Kernel Smoothing',
    author = 'Dennis Wylie',
    author_email = 'denniswylie@austin.utexas.edu',
    py_modules = ['sarks'],
    install_requires = ['biopython', 'editdistance', 'intervaltree', 'numpy', 'pandas', 'pyfaidx', 'scipy']
)
