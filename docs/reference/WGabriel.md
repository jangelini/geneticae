# Weighted GabrielEigen imputation method

agregar descripcion.

## Usage

``` r
WGabriel(DBmiss, Winf, Wsup)
```

## Arguments

- DBmiss:

  a data frame or matrix that contains the genotypes in the rows and the
  environments in the columns when there are no replications of the
  experiment.

- Winf:

  inferior weight

- Wsup:

  superior weight

## Value

A list containing:

- Peso: weight that provides the best predictive difference;

- NumeroIterWGabriel: the final number of iterations;

- CritConvergWGabriel: the maximum change of the estimated values for
  missing cells in the last step of iteration (the precision of
  convergence);

- GabrielWImput: the imputed matrix (filled in with the missing values
  estimated by the Weighted GabrielEigen procedure).

## References

Arciniegas-Alarcón S., García-Peña M., Krzanowski W.J., Dias C.T.S.
(2014). An alternative methodology for imputing missing data in trials
with genotype-byenvironment interaction: some new aspects. Biometrical
Letters 51, 75-88.
