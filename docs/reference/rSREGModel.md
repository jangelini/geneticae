# Site Regression model

The Site Regression model (also called genotype +
genotype-by-environment (GGE) model) is a powerful tool for effective
analysis and interpretation of data from multi-environment trials in
breeding programs. There are different functions in R to fit the SREG
model, however, this function has the following improvements:

- Includes recently published robust versions of the SREG model
  (Angelini et al., 2022).

- It can be used for data from trials with repetitions (there is no need
  to calculate means beforehand).

- Other variables not used in the analysis can be present in the
  dataset.

## Usage

``` r
rSREGModel(
  Data,
  genotype = "gen",
  environment = "env",
  response = "yield",
  rep = NULL,
  model = "SREG",
  SVP = "symmetrical"
)
```

## Arguments

- Data:

  dataframe with genotypes, environments, repetitions (if any) and the
  phenotypic trait of interest. Additional variables that will not be
  used in the model may be present in the data.

- genotype:

  column name for genotypes.

- environment:

  column name for environments.

- response:

  column name for the phenotypic trait.

- rep:

  column name for replications. If this argument is NULL, there are no
  replications in the data. Defaults to NULL.

- model:

  method for fitting the SREG model:
  \`"SREG"\`,\`"CovSREG"\`,\`"hSREG"\` or \`"ppSREG"\` (see References).
  Defaults to \`"SREG"\`.

- SVP:

  method for singular value partitioning. Either \`"row"\`,
  \`"column"\`, or \`"symmetrical"\`. Defaults to \`"symmetrical"\`.

## Value

A list of class `GGE_Model` containing:

- model:

  SREG model version.

- coordgenotype:

  plotting coordinates for each genotype in every component.

- coordenviroment:

  plotting coordinates for each environment in every component.

- eigenvalues:

  vector of eigenvalues for each component.

- vartotal:

  overall variance.

- varexpl:

  percentage of variance explained by each component.

- labelgen:

  genotype names.

- labelenv:

  environment names.

- axes:

  axis labels.

- Data:

  scaled and centered input data.

- SVP:

  name of SVP method.

A biplot of class `ggplot`

## Details

A linear model by robust regression using an M estimator proposed by
Huber (1964, 1973) fitted by iterated re-weighted least squares, in
combination with three robust SVD/PCA procedures, resulted in a total of
three robust SREG alternatives. The robust SVD/PCA considered were:

- CovSREG: robust PCA that is obtained by replacing the classical
  estimates of location and covariance by their robust analogues using
  Minimum Regularized Covariance Determinant (MRCD) approach;

- hSREG: robust PCA method that tries to combine the advantages of both
  approaches, PCA based on a robust covariance matrix and based on
  projection pursuit;

- ppSREG: robust PCA that uses the projection pursuit and directly
  calculates the robust estimates of the eigenvalues and eigenvectors
  without going through robust covariance estimation. It is a very
  attractive method for bigdata situations, which are very common in
  METs (a few genotypes tested in a large number of environments), as
  the principal components can be calculated sequentially.

## References

Julia Angelini, Gabriela Faviere, Eugenia Bortolotto, Gerardo Domingo
Lucio Cervigni & Marta Beatriz Quaglino (2022) Handling outliers in
multi-environment trial data analysis: in the direction of robust SREG
model, Journal of Crop Improvement,
<https://doi.org/10.1080/15427528.2022.2051217>

## Examples

``` r
 library(geneticae)

 # Data without replication
 library(agridat)
 data(yan.winterwheat)
 GGE1 <- rSREGModel(yan.winterwheat, genotype="gen", environment="env", response="yield")

 # Data with replication
 data(plrv)
 GGE2 <- rSREGModel(plrv, genotype = "Genotype", environment = "Locality",
                  response = "Yield", rep = "Rep")
```
