# Robust AMMI Model

Fits a classical or robust Additive Main effects and Multiplicative
Interaction (AMMI) model.

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

  a dataframe with genotypes, environments and the phenotypic trait.

- genotype:

  column name containing genotypes.

- environment:

  column name containing environments.

- response:

  column name containing the phenotypic trait.

- rep:

  column name containing replications. If provided, means are
  calculated.

- Ncomp:

  number of principal components to retain.

- type:

  method for fitting: \`"AMMI"\`, \`"rAMMI"\`, \`"hAMMI"\`, \`"gAMMI"\`,
  \`"lAMMI"\` or \`"ppAMMI"\`.
