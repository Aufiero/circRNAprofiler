
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![HitCount](http://hits.dwyl.io/Aufiero/circRNAprofiler.svg)](http://hits.dwyl.io/Aufiero/circRNAprofiler)

circRNAprofiler: an R-based computational framework for the downstream analysis of circular RNAs
================================================================================================

circRNAprofiler is an R-based framework that only requires an R installation and offers 15 modules for a comprehensive in silico analysis of circRNAs. This computational framework allows to combine and analyze circRNAs previously detected by multiple publicly available annotation-based circRNA detection tools. It covers different aspects of circRNAs analysis from differential expression analysis, evolutionary conservation, biogenesis to functional analysis. The pipeline used by circRNAprofiler is highly automated and customizable. Furthermore, circRNAprofiler includes additional functions for data visualization which facilitate the interpretation of the results.

Installation
------------

You can install circRNAprofiler from Bioconductor using:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("circRNAprofiler")
```

To download the development version on github use:

``` r
BiocManager::install("Aufiero/circRNAprofiler")
```

Documentation
-------------

Detailed information on using this package can be found in the package vignettes:

``` r
browseVignettes("circRNAprofiler")
```

Bugs and Issues
---------------

We work hard to ensure that circRNAprofiler is a powerful tool empowering your research. However, no software is free of bugs and issues, therefore we would love to get feedback from our users.
