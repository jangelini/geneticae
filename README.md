
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

Understanding the relationship between crops performance and environment
is a key problem for plant breeders and geneticists. In advanced stages
of breeding programs, where few genotypes are evaluated,
multi-environmental trials (MET) is one of the most common experiments.
They are conducted by testing a number of genotypes across multiple
environments, allowing the identification of superior genotypes. Crop
performance, the observed phenotype, is a function of genotype (G),
environment (E) and genotype x environment interaction (GEI). METs are
essential due to the presence of GEI which generates differential
genotypic responses in the different environments evaluated (Crossa et
al., 1990; Cruz Medina, 1992; Kang and Magari, 1996). This is a
particular problem for plant breeders (Giauffret et al., 2000),
therefore appropriate statistical methods should be used to obtain an
adequate GEI analysis.

geneticae package provides tools to analize data from advanced stages of
breeding programs, where few genotypes are evaluated. Among the
functions available in the package are the AMMI model and the SREG
model, and the biplots that are obtained from them. In addition,
functions are included that allow to overcome the fragility of the
methods between atypical observations, robust AMMI and imputation
techniques since the methods do not work in the presence of atypical
observations. Unlike the other existing functions that allow adjusting
these models, it is less restrictive in terms of the input data set.
Also, the biplots created by the package are ggplot objects, which means
that after a plot is created it can be further customized using various
functions from the ggplot2 package.

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
# install.packages('devtools')
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

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />