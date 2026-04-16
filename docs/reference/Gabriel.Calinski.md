# GabrielEigen imputation method

GabrielEigen imputation method

## Usage

``` r
Gabriel.Calinski(X)
```

## Arguments

- X:

  a data frame or matrix with genotypes in rows and environments in
  columns when there are no replications of the experiment.

## Value

A list containing:

- NumeroIterGabriel: final number of iterations;

- CritConvergGabriel: maximum change of the estimated values for missing
  cells in the last step of iteration (precision of convergence);

- GabrielImput: imputed matrix (filled in with missing values estimated
  by the GabrielEigein procedure).

## References

Arciniegas-Alarcón S., García-Peña M., Dias C.T.S., Krzanowski W.J.
(2010). *An alternative methodology for imputing missing data in trials
with genotype-by-environment interaction*. Biometrical Letters 47, 1–14.
