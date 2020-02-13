# sarks installation

## Requirements

- for R version: Java >= 1.8, R >= 3.3.0 with rJava successfully installed
  - relies on CRAN package
    [**remotes**](https://cran.r-project.org/web/packages/remotes/index.html)
    for installation
  - depends on CRAN packages
    [**rJava**](https://cran.r-project.org/web/packages/rJava/index.html),
    [**utils**](https://www.rdocumentation.org/packages/utils/versions/3.6.2),
    [**cluster**](https://cran.r-project.org/web/packages/cluster/index.html),
    [**binom**](https://cran.r-project.org/web/packages/binom/index.html)
  - and Bioconductor package
    [**IRanges**](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
- for Java version: Java >= 1.8

## Installation: R

1. Run the following R code within an R session:
   ```R
   ## if you don't already have remotes installed, uncomment and run:
   # install.packages('remotes')
   library(remotes)
   install_github('denniscwylie/sarks')
   ## alternatively, to build vignette as well, try uncommenting and running:
   # install_github('denniscwylie/sarks', build_vignettes=TRUE)
   
   ```

2. Test the installation by going through the simulated data example
   as described in the vignette (either the pdf version, if you built it
   while installing, or the [abridged markdown vignette](sarks_vignette.md)).

## Installation: Java

1. Copy sarks.jar from inst/java/ subdirectory of this repository
   to convenient location
   
2. Test the installation by going through the simulated data example
   as describedin [README.md](README.md)
   
