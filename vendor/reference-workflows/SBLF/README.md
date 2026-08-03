
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `SBLF`

[![](https://img.shields.io/badge/devel%20version-1.0.0.9000-blue.svg)](https://github.com/umich-biostatistics/SBLF)
[![](https://img.shields.io/github/languages/code-size/umich-biostatistics/SBLF.svg)](https://github.com/umich-biostatistics/SBLF)

## Overview

This `R` package provides a user-friendly interface for fitting spatial
Bayesian latent factor models. Uncover the latent brain image structure
using image on image regression.

## Installation

The package is currently only runnable on Linux. Windows or Mac users
can run a virtual Linux distribution on their machine or use a cluster.

To install, open `R` and run:

If the devtools package is not yet installed, install it first:

``` r
install.packages('devtools')
library(devtools)
```

``` r
# install SBLF from Github:
install_github('umich-biostatistics/SBLF') 
library(SBLF)
```

## Example usage

### Quick example:

``` r
?SBLF
mod = SBLF(xtrain, xtest, ztrain, ztest, voxel_loc, seed = 1234, burnin = 250, iter = 500)
```
