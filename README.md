# sarks

__SArKS__ (Suffix Array Kernel Smoothing) is an algorithm for
identifying sequence motifs correlated with numeric scores (such as
differential expression statistics from RNA-seq experiments). The
paper describing the algorithm may be found at:

https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797

A preprint of the article is also available on biorxiv at:

https://www.biorxiv.org/content/early/2018/10/25/133934


## Installation

SArKS is implemented in Java (1.8 or greater) with interactive use
facilitated through an R package built using
[**rJava**](https://cran.r-project.org/web/packages/rJava/index.html).

Once these dependencies have been installed and correctly configured,
you can install \R{sarks} by running the following code within an R
session:
```R
## if you don't already have remotes installed, uncomment and run:
# install.packages('remotes')
library(remotes)
install_github('denniscwylie/sarks')
## alternatively, to build vignette as well, try uncommenting and running:
# install_github('denniscwylie/sarks', build_vignettes=TRUE)
```

### Alternative installation: Java only

1. Copy sarks.jar from inst/java/ subdirectory of this repository
   to convenient location

2. Test the installation by going through the simulated data example
   using sarks.jar as described below

## Using sarks

This project implements the SArKS algorithm in the java package
contained in sarks.jar, which can also be run as part of the R package
sarks.

### Using the R package sarks

For most users, we would recommend trying out the R package, which can
be installed as described above.

The sarks vignette is the best place to start to learn how to use the
R version of sarks.

The full vignette is available as a pdf if you use the
`"build_vignettes=TRUE"` option when installing sarks in R; otherwise,
you can take a look at the [abridged markdown vignette](sarks_vignette.md).

### Direct command-line usage of jar file

For detailed information on command-line usage of sarks.jar and
associated scripts, consult [user_guide.md](user_guide.md).

The best way to learn how to use sarks is to read through the example
scripts

```
examples/*_example.sh
```

(markdown versions of each of the examples are available as well)
included in the github repository.

These examples are taken from the data sets analyzed in the SArKS
paper, including the toy simulated data set as well as the analyses of
the upstream (5' of transcription start site) and downstream (3' of
transcription start site) DNA regions for mouse genes whose expression
profiles were quantified in the studies:

- Mo, Alisa, et al. "Epigenomic signatures of neuronal diversity in the
  mammalian brain." Neuron 86.6 (2015): 1369-1384.
- Close, Jennie L., et al. "Single-cell profiling of an in vitro model
  of human interneuron development reveals temporal dynamics of cell
  type production and maturation." Neuron 93.5 (2017): 1035-1048.


### Simulated data example

The simulated data set consists of the 30 sequences contained in

- examples/simulated_seqs.fa

together with the associated scores contained in

- examples/simulated_scores.tsv

The file

[examples/simulated_example.md](examples/simulated_example.md)

uses the utility scripts also contained in the examples folder to
analyze these sequences and scores. After moving to the examples
directory,

```
cd examples/
```

I recommend reading through the example and running the commands
contained within individually at the command line as you get to them.


### Mo 2015 downstream example

After going through the simulated example, try sarks out on the Mo
2015 downstream seqs. An example of how to do this can be found in the

[examples/mo2015\_downstream\_example.md](examples/mo2015_downstream_example.md)

file; again I would recommend reading through the example and running
the commands line-by-line as you get to them.

### Mo 2015 upstream example

The [Mo 2015 upstream example](examples/mo2015_upstream_example.md) is
similar to the downstream example but
1. runs on a larger data set and
2. uses a few more SArKS features, including spatial smoothing.

Because of these features, this example will run slower and require
a good deal more memory than the other examples included.
