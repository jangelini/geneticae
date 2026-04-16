# AMMI Biplots with ggplot2

Produces a biplot for objects of class 'AMMI'.

## Usage

``` r
rAMMIPlot(
  model_res,
  colGen = "gray47",
  colEnv = "darkred",
  sizeGen = 6,
  sizeEnv = 6,
  titles = TRUE,
  footnote = TRUE,
  axis_expand = 1.2,
  limits = TRUE,
  axes = TRUE,
  axislabels = TRUE
)
```

## Arguments

- model_res:

  an object of class 'AMMI' from AMMIModel.

- colGen:

  genotype colour. Defaults to "gray47".

- colEnv:

  environment colour. Defaults to "darkred".

- sizeGen:

  genotype text size.

- sizeEnv:

  environment text size.

- titles:

  logical, show plot title.

- footnote:

  logical, show footnote with explained variance.

- axis_expand:

  expansion factor for axis limits.

- limits:

  logical. If \`TRUE\` axes are automatically rescaled. Defaults to
  \`TRUE\`.

- axes:

  logical, if this argument is \`TRUE\` axes passing through the origin
  are drawn. Defaults to \`TRUE\`.

- axislabels:

  logical, if this argument is \`TRUE\` labels axes are included.
  Defaults to \`TRUE\`
