# sarks

__SArKS__ (Suffix Array Kernel Smoothing) is an algorithm for
identifying sequence motifs correlated with numeric scores (such as
differential expression statistics from RNA-seq experiments). A
preprint describing the algorithm may be found at:

https://www.biorxiv.org/content/early/2018/03/11/133934


## Installation
See [INSTALL.md](INSTALL.md).


## Using sarks

This project implements the SArKS algorithm in the python module
sarks. To facilitate usage we have also provided several utility
shell scripts making use of this python module, including:

- examples/sarkselect.py for identifying k-mer motifs using SArKS
- examples/sarkstest.py for estimating false positive rates
- examples/extract_kmers.py for quickly pulling list of sarks k-mers
  from sarkselect.py output files
- examples/cluster_seqs.R for (optionally) clustering sarks k-mers
- examples/count_kmers.py for counting or locating sarks k-mers
  (or clusters of k-mers) in sequences contained in a fasta file

The best way to learn how to use sarks is to read through the example
scripts provided in the examples/ folder included in the github
repository.

These examples are taken from the data sets analyzed in the SArKS
paper, including the toy simulated data set as well as the analyses of
the upstream (5' of transcription start site) and downstream (3' of
transcription start site) DNA regions for mouse genes whose expression
profiles were quantified in the studies:

Mo, Alisa, et al. "Epigenomic signatures of neuronal diversity in the
mammalian brain." Neuron 86.6 (2015): 1369-1384.

Close, Jennie L., et al. "Single-cell profiling of an in vitro model
of human interneuron development reveals temporal dynamics of cell
type production and maturation." Neuron 93.5 (2017): 1035-1048.


### Simulated data example

The simulated data set consists of the 30 sequences contained in

examples/simulated_seqs.fa

together with the associated scores contained in

examples/simulated_scores.tsv

The shell script

examples/simulated_example.sh

uses the small utility scripts also contained in the examples folder
to analyze these sequences and scores. After moving to the examples
directory,
```
cd examples/
```
I recommend reading through the script file and running the commands
contained within individually at the shell script as you get to them.


### Mo 2015 downstream example

After going through the simulated example, try sarks out on the Mo
2015 downstream seqs. Recommend following the same procedure as for
the simulated data set but referencing the mo2015 downstream example
shell script file; also note that you will have to gunzip the fasta
file first:

```
gunzip mo2015_downstream_seqs.fa.gz
```

NOTE: as this is a much larger data set than the simulated set,
running this example will take a bit longer. Only one smoothing window
(and no spatial smoothing) is applied in this analysis so that motif
discovery can be run in a minute or two. The number of permutations
used in initial permutational threshold setting has also been reduced
in the interest of speed, so results will vary a bit more from run to
run.

The sarkstest.py step applied to estimate false positive rates is
still set to use 250 permutations and will thus likely take > 10
minutes to complete depending on hardware.


### Mo 2015 upstream example

NOTE: UPSTREAM EXAMPLE REQUIRES > 32 GIGABYTES AVAILABLE RAM

(recommend at least 48G available to run comfortably)

The Mo 2015 upstream example script is similar to the downstream
example but (1)) runs on a larger data set and (2) uses a few more
SArKS features, including spatial smoothing.

Because of these features, this example will run slower and require
a good deal more memory than the other examples included.
