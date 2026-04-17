# Imputation of missing cells in two-way data sets

Missing values are not allowed by the AMMI, GGE or SREG methods. This
function provides several methods to impute missing observations in data
from multi-environment trials and to subsequently adjust the mentioned
methods.

## Usage

``` r
imputation(
  Data,
  genotype = "gen",
  environment = "env",
  response = "yield",
  rep = NULL,
  type = "EM-AMMI",
  nPC = 2,
  initial.values = NA,
  precision = 0.01,
  maxiter = 1000,
  change.factor = 1,
  simplified.model = FALSE,
  scale = TRUE,
  method = "EM",
  row.w = NULL,
  coeff.ridge = 1,
  seed = NULL,
  nb.init = 1,
  Winf = 0.8,
  Wsup = 1
)
```

## Arguments

- Data:

  dataframe containing genotypes, environments, repetitions (if any) and
  the phenotypic trait of interest. Other variables that will not be
  used in the analysis can be present.

- genotype:

  column name containing genotypes.

- environment:

  column name containing environments.

- response:

  column name containing the phenotypic trait.

- rep:

  column name containing replications. If this argument is NULL, there
  are no replications available in the data. Defaults to NULL.

- type:

  imputation method. Either "EM-AMMI", "EM-GGE", "EM-SREG", "EM-bSREG",
  "Gabriel", "Eigenvector", "WGabriel", "EM-PCA". Defaults to "EM-AMMI".

- nPC:

  number of components used to predict the missing values. Default to 2.

- initial.values:

  initial values of the missing cells. It can be a single value or a
  vector of length equal to the number of missing cells.

- precision:

  threshold for assessing convergence.

- maxiter:

  maximum number of iteration for the algorithm.

- change.factor:

  When \`change.factor\` is equal to 1, the previous approximation is
  changed with the new values (standard EM). Smaller values can help
  convergence if changes are cyclic.

- simplified.model:

  logical. If TRUE, calculates effects only in the first iteration to
  speed up convergence or help in cases where the regular procedure
  fails.

- scale:

  boolean. By default TRUE for "EM-PCA".

- method:

  "Regularized" or "EM" for "EM-PCA".

- row.w:

  row weights for "EM-PCA".

- coeff.ridge:

  ridge coefficient for "EM-PCA".

- seed:

  integer for random initialization in "EM-PCA".

- nb.init:

  number of random initializations for "EM-PCA".

- Winf:

  lower weight for WGabriel.

- Wsup:

  upper weight for WGabriel.

## Value

A matrix of the imputed data.

## Details

Often, multi-environment experiments are unbalanced because several
genotypes are not tested in some environments. Several methodologies
have been proposed in order to solve this lack of balance caused by
missing values, some of which are included in this function:

- EM-AMMI: an iterative scheme built round the above procedure is used
  to obtain AMMI imputations from the EM algorithm. The additive
  parameters are initially set by computing the grand mean, genotype
  means and environment means obtained from the observed data. The
  residuals for the observed cells are initialized as the cell mean
  minus the genotype mean minus the environment mean plus the grand
  mean, and interactions for the missing positions are initially set to
  zero. The initial multiplicative parameters are obtained from the SVD
  of this matrix of residuals, and the missing values are filled by the
  appropriate AMMI estimates. In subsequent iterations, the usual AMMI
  procedure is applied to the completed matrix and the missing values
  are updated by the corresponding AMMI estimates. The arguments used
  for this method are:initial.values, precision, maxiter, change.factor
  and simplified.model

- EM-GGE: Iterative SVD-based imputation focusing on G+GE.

- EM-SREG: Iterative algorithm using the Sites Regression model.
  Supports variants like standard SVD and Bayesian PCA (EM-bSREG).

- Gabriel: combines regression and lower-rank approximation using SVD.
  This method initially replaces the missing cells by arbitrary values,
  and subsequently the imputations are refined through an iterative
  scheme that defines a partition of the matrix for each missing value
  in turn and uses a linear regression of columns (or rows) to obtain
  the new imputation. The arguments used for this method is only the
  dataframe.

- WGabriel: is a a modification of Gabriel method that uses weights
  chosen by cross-validation. The arguments used for this method are
  Winf and Wsup.

- EM-PCA: impute the missing entries of a mixed data using the iterative
  PCA algorithm. The algorithm first consists imputing missing values
  with initial values. The second step of the iterative PCA algorithm is
  to perform PCA on the completed dataset to estimate the parameters.
  Then, it imputes the missing values with the reconstruction formulae
  of order nPC (the fitted matrix computed with nPC components for the
  scores and loadings). These steps of estimation of the parameters via
  PCA and imputation of the missing values using the fitted matrix are
  iterate until convergence. The arguments used for this methods are:
  nPC, scale, method, row.w, coeff.ridge, precision, seed, nb.init and
  maxiter

## References

Paderewski, J. (2013). *An R function for imputation of missing cells in
two-way data sets by EM-AMMI algorithm*. Communications in Biometry and
Crop Science 8, 60–69.

Yan, W. (2013). *Biplot analysis of incomplete two-way data*. Crop
Science, 53(1), 48-57.
[doi:10.2135/cropsci2012.05.0301](https://doi.org/10.2135/cropsci2012.05.0301)

Arciniegas-Alarcón, S., García-Peña, M., Krzanowski, W., & Dias, C. T.
S. (2014b). *An alternative methodology for imputing missing data in
trials with genotype-by-environment interaction: some new aspects*.
Biometrical Letters, 51(2), 75-88.
[doi:10.2478/bile-2014-0006](https://doi.org/10.2478/bile-2014-0006)

Angelini, J., Cervigni, G. D. L., & Quaglino, M. B. (2024). *New
imputation methodologies for genotype-by-environment data: an extensive
study of properties of estimators*. Euphytica, 220(6), 92.
[doi:10.1007/s10681-024-03344-z](https://doi.org/10.1007/s10681-024-03344-z)

Julie Josse, Francois Husson (2016). missMDA: A Package for Handling
Missing Values in Multivariate Data Analysis. Journal of Statistical
Software 70, 1-31.

Arciniegas-Alarcón S., García-Peña M., Dias C.T.S., Krzanowski W.J.
(2010). *An alternative methodology for imputing missing data in trials
with genotype-by-environment interaction*. Biometrical Letters 47, 1–14.

Arciniegas-Alarcón S., García-Peña M., Krzanowski W.J., Dias C.T.S.
(2014). *An alternative methodology for imputing missing data in trials
with genotype-byenvironment interaction: some new aspects.* Biometrical
Letters 51, 75-88.

## Examples

``` r
library(geneticae)
# Data without replications
library(agridat)
data(yan.winterwheat)

# generating missing values
yan.winterwheat[1,3]<-NA
yan.winterwheat[3,3]<-NA
yan.winterwheat[2,3]<-NA

imputation(yan.winterwheat, genotype = "gen", environment = "env",
           response = "yield", type = "EM-AMMI")
#>         BH93  EA93  HW93  ID93  KE93  NN93  OA93  RN93  WP93
#> Ann 4.150120 4.150 2.849 3.084 5.940 4.450 4.351 4.039 2.672
#> Ari 4.035814 4.771 2.912 3.506 5.699 5.152 4.956 4.386 2.938
#> Aug 4.305244 4.578 3.098 3.460 6.070 5.025 4.730 3.900 2.621
#> Cas 4.732000 4.745 3.375 3.904 6.224 5.340 4.226 4.893 3.451
#> Del 4.390000 4.603 3.511 3.848 5.773 5.421 5.147 4.098 2.832
#> Dia 5.178000 4.475 2.990 3.774 6.583 5.045 3.985 4.271 2.776
#> Ena 3.375000 4.175 2.741 3.157 5.342 4.267 4.162 4.063 2.032
#> Fun 4.852000 4.664 4.425 3.952 5.536 5.832 4.168 5.060 3.574
#> Ham 5.038000 4.741 3.508 3.437 5.960 4.859 4.977 4.514 2.859
#> Har 5.195000 4.662 3.596 3.759 5.937 5.345 3.895 4.450 3.300
#> Kar 4.293000 4.530 2.760 3.422 6.142 5.250 4.856 4.137 3.149
#> Kat 3.151000 3.040 2.388 2.350 4.229 4.257 3.384 4.071 2.103
#> Luc 4.104000 3.878 2.302 3.718 4.555 5.149 2.596 4.956 2.886
#> m12 3.340000 3.854 2.419 2.783 4.629 5.090 3.281 3.918 2.561
#> Reb 4.375000 4.701 3.655 3.592 6.189 5.141 3.933 4.208 2.925
#> Ron 4.940000 4.698 2.950 3.898 6.063 5.326 4.302 4.299 3.031
#> Rub 3.786000 4.969 3.379 3.353 4.774 5.304 4.322 4.858 3.382
#> Zav 4.238000 4.654 3.607 3.914 6.641 4.830 5.014 4.363 3.111

# Data with replications
data(plrv)
head(plrv)
#>   Genotype Locality Rep WeightPlant WeightPlot    Yield
#> 1   102.18     Ayac   1   0.5100000       5.10 18.88889
#> 2   104.22     Ayac   1   0.3450000       2.76 12.77778
#> 3   121.31     Ayac   1   0.5425000       4.34 20.09259
#> 4   141.28     Ayac   1   0.9888889       8.90 36.62551
#> 5   157.26     Ayac   1   0.6250000       5.00 23.14815
#> 6    163.9     Ayac   1   0.5120000       2.56 18.96296
plrv$Yield[plrv$Locality == "Ayac" & plrv$Rep %in% c(1, 2, 3) & plrv$Genotype == '102.18'] <- NA

imputation(plrv, nPC = 2,genotype = "Genotype", environment = "Locality", 
           response = "Yield", rep ='Rep', type = "EM-AMMI")
#>             Ayac    Hyo-02    LM-02    LM-03     SR-02     SR-03
#> 102.18  23.55025 28.888889 32.03704 46.77778 13.518519 11.769547
#> 104.22  21.45062 53.518519 39.19753 50.41838 16.049383  7.098765
#> 121.31  23.46022 41.296296 38.39506 63.70370  2.500000 11.255144
#> 141.28  31.84401 60.462963 33.95062 77.56790 19.234568 15.477366
#> 157.26  19.66980 41.388889 45.16049 76.98553 23.950617 14.555556
#> 163.9   17.53792 29.537037 28.88889 32.02949 12.716049  7.795414
#> 221.19  15.41358 32.037037 31.02469 43.45267  8.543210  7.437586
#> 233.11  24.28326 50.555556 29.19753 47.33333 13.148148  7.481481
#> 235.6   29.91358 73.518519 40.37037 56.13580 14.802469 17.067901
#> 241.2   20.44444 36.018519 35.74074 46.24691 11.086420  8.505291
#> 255.7   26.06702 47.037037 32.67901 40.79244 19.185185 17.777778
#> 314.12  17.32510 49.444444 34.50617 57.24588  8.530864  1.987654
#> 317.6   26.61376 53.425926 42.34568 64.01235 14.814815 10.742455
#> 319.20  25.77503 56.666667 32.96296 86.80802 21.209877  9.123457
#> 320.16  30.32922 31.111111 35.28395 43.29012 13.629630  4.444444
#> 342.15  19.89712 40.740741 27.53086 38.80033 17.592593 11.518519
#> 346.2   21.57476 32.685185 25.55556 32.03704 18.024691 13.173280
#> 351.26  31.74897 50.185185 29.25926 72.02407 20.370370 13.106996
#> 364.21  26.63933 52.407407 37.90123 57.06584 13.518519 16.826132
#> 402.7   19.29698 42.500000 31.23457 49.76543 12.839506  9.228395
#> 405.2   28.66735 35.740741 32.34568 43.25926 16.790123 17.116598
#> 406.12  19.58652 59.814815 37.77778 53.58848 13.827160 11.504605
#> 427.7   26.08907 55.648148 44.44444 58.33608 21.234568 11.388889
#> 450.3   28.72428 50.185185 36.88889 72.24198 15.432099 13.703704
#> 506.2   25.00000 46.759259 45.55556 53.24966 18.148148 10.884774
#> Canchan 21.32716 47.777778 21.60494 59.24691  9.629630  2.421125
#> Desiree 18.76543  8.888889 20.37037 27.42747 10.061728 11.420243
#> Unica   21.30144 72.222222 47.83951 57.53519 18.246914 17.478738
           
imputation(plrv, genotype = "Genotype", environment = "Locality", 
           response = "Yield", rep ='Rep', type = "EM-SREG")
#>             Ayac    Hyo-02    LM-02    LM-03     SR-02     SR-03
#> 102.18  21.67165 28.888889 32.03704 46.77778 13.518519 11.769547
#> 104.22  21.45062 53.518519 39.19753 50.41838 16.049383  7.098765
#> 121.31  23.46022 41.296296 38.39506 63.70370  2.500000 11.255144
#> 141.28  31.84401 60.462963 33.95062 77.56790 19.234568 15.477366
#> 157.26  19.66980 41.388889 45.16049 76.98553 23.950617 14.555556
#> 163.9   17.53792 29.537037 28.88889 32.02949 12.716049  7.795414
#> 221.19  15.41358 32.037037 31.02469 43.45267  8.543210  7.437586
#> 233.11  24.28326 50.555556 29.19753 47.33333 13.148148  7.481481
#> 235.6   29.91358 73.518519 40.37037 56.13580 14.802469 17.067901
#> 241.2   20.44444 36.018519 35.74074 46.24691 11.086420  8.505291
#> 255.7   26.06702 47.037037 32.67901 40.79244 19.185185 17.777778
#> 314.12  17.32510 49.444444 34.50617 57.24588  8.530864  1.987654
#> 317.6   26.61376 53.425926 42.34568 64.01235 14.814815 10.742455
#> 319.20  25.77503 56.666667 32.96296 86.80802 21.209877  9.123457
#> 320.16  30.32922 31.111111 35.28395 43.29012 13.629630  4.444444
#> 342.15  19.89712 40.740741 27.53086 38.80033 17.592593 11.518519
#> 346.2   21.57476 32.685185 25.55556 32.03704 18.024691 13.173280
#> 351.26  31.74897 50.185185 29.25926 72.02407 20.370370 13.106996
#> 364.21  26.63933 52.407407 37.90123 57.06584 13.518519 16.826132
#> 402.7   19.29698 42.500000 31.23457 49.76543 12.839506  9.228395
#> 405.2   28.66735 35.740741 32.34568 43.25926 16.790123 17.116598
#> 406.12  19.58652 59.814815 37.77778 53.58848 13.827160 11.504605
#> 427.7   26.08907 55.648148 44.44444 58.33608 21.234568 11.388889
#> 450.3   28.72428 50.185185 36.88889 72.24198 15.432099 13.703704
#> 506.2   25.00000 46.759259 45.55556 53.24966 18.148148 10.884774
#> Canchan 21.32716 47.777778 21.60494 59.24691  9.629630  2.421125
#> Desiree 18.76543  8.888889 20.37037 27.42747 10.061728 11.420243
#> Unica   21.30144 72.222222 47.83951 57.53519 18.246914 17.478738
```
