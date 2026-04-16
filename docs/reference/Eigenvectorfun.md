# Eigenvector Imputation Function (Internal)

Internal function for GxE imputation using the Krzanowski (1988)
eigenvector approach with a leave-one-out strategy.

## Usage

``` r
Eigenvectorfun(X, f)
```

## Arguments

- X:

  A matrix with missing values (NAs).

- f:

  Number of components (rank) to use for the reconstruction.

## Value

A list with the number of iterations, convergence status, final rank
used, and the imputed matrix.
