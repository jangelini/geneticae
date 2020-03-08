
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geneticae

# <img src="man/figures/baseplot.png" align="right">

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/geneticae)](https://CRAN.R-project.org/package=geneticae)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Downloads](https://cranlogs.r-pkg.org/badges/geneticae?color=blue)](https://cran.rstudio.com/package=geneticae)
[![Codecov test
coverage](https://codecov.io/gh/r-lib/geneticae/branch/master/graphs/badge.svg)](https://codecov.io/gh/r-lib/geneticae?branch=master)
<!-- badges: end -->

The geneticae package provides tools to analize data from advanced
stages of breeding programs, where few genotypes are evaluated.

geneticae package there are several functions for use after fitting
models (AMMI or SREG). The plots created by the package are ggplot
objects, which means that after a plot is created it can be further
customized using various functions from the ggplot2 package.

## Getting Started

If you are just getting started with geneticae we recommend starting
with the tutorial
[*vignettes*](file:///F:/Especializacion%20en%20bioinformatica/Para%20proyecto%20final/Geneticae%20Package/geneticae/docs/articles/vignettes.html),
and the examples throughout the package
[*documentation*](file:///F:/Especializacion%20en%20bioinformatica/Para%20proyecto%20final/Geneticae%20Package/geneticae/docs/reference/index.html).

## Installation

You can install the released version of geneticae from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("geneticae")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jangelini/geneticae")
```

## Example

Some quick examples obtained from the geneticae package.

``` r
library(geneticae)

data(yan.winterwheat)
GGE1 <- GGEmodel(yan.winterwheat, centering = "tester")
GGEPlot(GGE1)
```
