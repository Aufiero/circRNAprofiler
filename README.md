
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FAufiero%2FcircRNAprofiler&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

# circRNAprofiler: an R-based computational framework for the downstream analysis of circular RNAs

circRNAprofiler is an R-based framework that only requires an R
installation and offers 15 modules for a comprehensive in silico
analysis of circRNAs. This computational framework allows to combine and
analyze circRNAs previously detected by multiple publicly available
annotation-based circRNA detection tools. It covers different aspects of
circRNAs analysis from differential expression analysis, evolutionary
conservation, biogenesis to functional analysis. The pipeline used by
circRNAprofiler is highly automated and customizable. Furthermore,
circRNAprofiler includes additional functions for data visualization
which facilitate the interpretation of the results.

## Installation

You can install the latest release of circRNAprofiler from Bioconductor
using:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("circRNAprofiler")
```

To install the latest development version use:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("circRNAprofiler")
```

## Documentation

Detailed information on using this package can be found in the package
vignettes:

``` r
browseVignettes("circRNAprofiler")
```

## Bugs and Issues

We work hard to ensure that circRNAprofiler is a powerful tool
empowering your research. However, no software is free of bugs and
issues, therefore we would love to get feedback from our users.

## Citation

If you find this code useful in your research, please cite:

Aufiero, S., Reckman, Y.J., Tijsen, A.J. et al. circRNAprofiler: an
R-based computational framework for the downstream analysis of circular
RNAs. BMC Bioinformatics 21, 164 (2020).
<https://doi.org/10.1186/s12859-020-3500-3>
