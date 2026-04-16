# Biplot Imputation Method

This function implements the Biplot imputation method as proposed by Yan
(2013). It is an iterative algorithm that uses the singular value
decomposition (SVD) to impute missing values in a genotype by
environment matrix.

## Usage

``` r
BiplotImputfun(X, precision = 0.01, max.iter = 1000, n_pc = 2)
```

## Arguments

- X:

  A data frame or matrix with genotypes in rows and environments in
  columns.

- precision:

  (optional) Convergence threshold. The algorithm stops when the
  relative change in imputed values is less than this value. Default is
  0.01.

- max.iter:

  (optional) Maximum number of iterations. Default is 1000.

- n_pc:

  Number of principal components to use for imputation. Default is 2.

## Value

A list containing:

- `X_imputed`: The final matrix with missing values filled.

- `iteration`: Total number of iterations performed.

- `convergence`: Final relative change reached.

- `fitted`: The fitted values from the AMMI/Biplot model.

## References

Yan, W. (2013). *Biplot analysis of incomplete two-way data*. Crop
Science, 53(1), 48-57.
[doi:10.2135/cropsci2012.05.0301](https://doi.org/10.2135/cropsci2012.05.0301)

Arciniegas-Alarcón, S., García-Peña, M., Krzanowski, W., & Dias, C. T.
S. (2014b). *An alternative methodology for imputing missing data in
trials with genotype-by-environment interaction: some new aspects*.
Biometrical Letters, 51(2), 75-88.
[doi:10.2478/bile-2014-0006](https://doi.org/10.2478/bile-2014-0006)
