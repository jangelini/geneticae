# EM-SREG Imputation Method

Iterative algorithm to impute missing values in two-way tables using the
Sites Regression (SREG) model. It supports several variants including
standard SVD and Bayesian PCA.

## Usage

``` r
EM.SREG(
  X,
  PC.nb = 1,
  initial.values = NA,
  precision = 0.01,
  max.iter = 1000,
  change.factor = 1,
  simplified.model = FALSE,
  type = c("EM-SREG", "EM-bSREG")
)
```

## Arguments

- X:

  A data frame or matrix with genotypes in rows and environments in
  columns.

- PC.nb:

  Number of principal components to be used. Default is 1.

- initial.values:

  (optional) Initial values for missing cells. If NA, initial values are
  obtained from column means (environment effects).

- precision:

  Convergence threshold. Default is 0.01.

- max.iter:

  Maximum number of iterations. Default is 1000.

- change.factor:

  Step size for updating missing values (standard is 1).

- simplified.model:

  Logical. If TRUE, effects are only calculated in the first iteration.

- type:

  Method type: "EM-SREG" (Standard), "EM-bSREG" (Bayesian).

## Value

A list containing:

- `X`: The imputed matrix.

- `iter`: The number of iterations until convergence.

## References

Angelini, J., Cervigni, G. D. L., & Quaglino, M. B. (2024). *New
imputation methodologies for genotype-by-environment data: an extensive
study of properties of estimators*. Euphytica, 220(6), 92.
[doi:10.1007/s10681-024-03344-z](https://doi.org/10.1007/s10681-024-03344-z)
