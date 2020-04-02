## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("jangelini/geneticae")

## ---- tidy=TRUE---------------------------------------------------------------
library(geneticae)
data(yan.winterwheat)
head(yan.winterwheat)

## ---- tidy=TRUE---------------------------------------------------------------
data(plrv)
head(plrv)

## ---- tidy=TRUE---------------------------------------------------------------
GGE1 <- GGEmodel(yan.winterwheat, genotype="gen", environment="env", response="yield")

## ---- tidy=TRUE---------------------------------------------------------------
names(GGE1)

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Biplot")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Selected Environment", selectedE = "OA93")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Selected Genotype", selectedG = "Kat")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Relationship Among Environments")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Comparison of Genotype", selectedG1 = "Kat", selectedG2 = "Cas")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Which Won Where/What")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Discrimination vs. representativeness")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Ranking Environments")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Ranking Genotypes")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE-------------
GGEPlot(GGE1, type="Mean vs. Stability")

