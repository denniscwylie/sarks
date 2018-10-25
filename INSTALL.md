# sarks installation

Requirements
------------
- python3.4 or later
  - setup.py should install required python modules if necessary
    (biopython, editdistance, intervaltree, numpy, pandas, pyfaidx, scipy)
- examples/cluster_seqs.R script uses philentropy and msa R libraries
- SeqAn C++ library (https://www.seqan.de/) version 1.4.x
  - on Ubuntu 16.04 or 18.04, can be installed using
    apt-get install seqan-dev
  - otherwise can be downloaded from
    packages.seqan.de/seqan-library/seqan-library-1.4.2.tar.bz2
  - if installed in a local directory, may need to edit the
    setup.py script to indicate location (see below)
- GNU make for compiling suffix-array.cpp and windginiimp.cpp utilities
- installation has been tested on unix systems, may need to be modified for others

Installation
------------

1. Clone the repository:
   ```bash
   git clone https://github.com/denniscwylie/sarks.git
   ```

2. Move to source directory:
   ```bash
   cd sarks
   ```

3. *If* SeqAn C++ library (v1.4) is installed in local directory:
   Open setup.py and edit the line
   incdir = ''
   replacing the empty string '' with the correct location on your system
   (a commented-out example of how to do so is shown in setup.py).

4. Install using pip3 to run setup.py; recommend building local and editable version
   ```bash
   pip3 install --user -e .
   ```
   If all goes well this should both set up the python module sarks and
   run make to build the utilities suffix-array and windginiimp.
   
   These two executables should have shown up in the directory to which the python code
   ```python
   site.USER_BASE + '/bin'
   ```
   evaluates (likely `$HOME/.local/bin`).

5. Make sure that the path for the directory indicated by the python code
   ```python
   site.USER_BASE + '/bin'
   ```
   is in your PATH. This is necessary so that the suffix-array and windginiimp
   utilities can be run by sarks.

6. Test the installation by going through the simulated data example
   as described in README.md.
