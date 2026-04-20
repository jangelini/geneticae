# Robust AMMI Model Fitting

Fits a classical or robust Additive Main effects and Multiplicative
Interaction (AMMI) model for genotype-by-environment data.

## Usage

``` r
rAMMIModel(
  Data,
  genotype = "gen",
  environment = "env",
  response = "Y",
  rep = NULL,
  Ncomp = 2,
  type = "AMMI"
)
```

## Arguments

- Data:

  a dataframe with genotypes, environments, repetitions (if any) and the
  phenotypic trait of interest. Other variables that will not be used in
  the analysis can be included.

- genotype:

  column name containing genotypes. Defaults to \`"gen"\`.

- environment:

  column name containing environments. Defaults to \`"env"\`.

- response:

  column name containing the phenotypic trait of interest. Defaults to
  \`"Y"\`.

- rep:

  column name containing replications. If this argument is \`NULL\`
  (default), it is assumed that the data already contains means per
  genotype in each environment. If provided, means are calculated
  automatically.

- Ncomp:

  number of principal components to retain for the interaction part.
  Defaults to 2.

- type:

  method for fitting the AMMI model: \`"AMMI"\` (classical),
  \`"rAMMI"\`, \`"hAMMI"\`, \`"gAMMI"\`, \`"lAMMI"\` or \`"ppAMMI"\`
  (robust variants). Defaults to \`"AMMI"\`.

## Value

A list of class `rAMMI` containing:

- gen_scores:

  Matrix of genotype scores (U \* D).

- env_scores:

  Matrix of environment loadings (V).

- eigenvalues:

  Vector of singular values for the retained components.

- gen_labels:

  Names of the genotypes.

- env_labels:

  Names of the environments.

- Ncomp:

  Number of principal components used.

- type:

  The fitting method used.

- vartotal:

  Total variance explained by the multiplicative terms.

## Details

To overcome the problem of data contamination with outlying
observations, Rodrigues, Monteiro and Lourenco (2015) propose a robust
AMMI model based on the M-Huber estimator and robust SVD/PCA procedures.

The \`type\` argument allows choosing between several robust strategies:

- **AMMI**: Classical AMMI model using Least Squares and standard SVD.

- **rAMMI**: Uses the L1 norm instead of the L2 norm to compute a robust
  approximation to the SVD (via pcaMethods).

- **hAMMI**: Uses the Hubert's approach (PcaHubert) combining
  projection-pursuit and robust covariance estimation.

- **gAMMI**: Uses the Grid search algorithm for PCA (PcaGrid).

- **lAMMI**: Performs PCA on the data projected onto a unit sphere
  (PcaLocantore).

- **ppAMMI**: Uses projection-pursuit (PcaProj) to calculate robust
  eigenvalues and eigenvectors.

## References

Rodrigues P.C., Monteiro A., Lourenco V.M. (2015). A robust AMMI model
for the analysis of genotype-by-environment data. Bioinformatics 32,
58-66.

## Examples

``` r
library(agridat)
data(yan.winterwheat)

# Classical AMMI
mod_ammi <- rAMMIModel(yan.winterwheat, genotype = "gen",
                       environment = "env", response = "yield", type = "AMMI")

# Robust AMMI (using Hubert's method)
mod_rammi <- rAMMIModel(yan.winterwheat, genotype = "gen",
                        environment = "env", response = "yield", type = "hAMMI")

```
