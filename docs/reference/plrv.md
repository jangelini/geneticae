# Clones from the PLRV population

resistance study to PLRV (Patato Leaf Roll Virus) causing leaf curl. 28
genotypes were experimented at 6 locations in Peru. Each clone was
evaluated three times in each environment, and yield, plant weight and
plot were registered.

## Usage

``` r
data(plrv)
```

## Format

Data frame with 504 observations and 6 variables (genotype, locality,
repetition, weightPlant, weightPlot and yield).

## References

Felipe de Mendiburu (2020). agricolae: Statistical Procedures for
Agricultural Research. R package version 1.3-2.
<https://CRAN.R-project.org/package=agricolae>

## Examples

``` r
library(geneticae)
data(plrv)
str(plrv)
#> 'data.frame':    504 obs. of  6 variables:
#>  $ Genotype   : Factor w/ 28 levels "102.18","104.22",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ Locality   : Factor w/ 6 levels "Ayac","Hyo-02",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Rep        : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ WeightPlant: num  0.51 0.345 0.542 0.989 0.625 ...
#>  $ WeightPlot : num  5.1 2.76 4.34 8.9 5 2.56 2.48 10.1 8.25 4.88 ...
#>  $ Yield      : num  18.9 12.8 20.1 36.6 23.1 ...
```
