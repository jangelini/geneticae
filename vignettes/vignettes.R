## ----eval=FALSE, tidy=TRUE-----------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("jangelini/geneticae")

## ---- tidy=TRUE----------------------------------------------------------
library(geneticae)
data(yan.winterwheat)
dat_yan <- yan.winterwheat
head(dat_yan)

## ---- tidy=TRUE----------------------------------------------------------
data(plrv)
dat_rep <- plrv
head(dat_rep)

## ---- tidy=TRUE----------------------------------------------------------
GGE1 <- GGEmodel(dat_yan, genotype="gen",environment="env", response="yield", centering = "tester")

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  GGE1_rep <- GGEmodel(dat_rep, genotype="Genotype",environment="Locality", response="Yield", rep="Rep", centering = "tester")

## ---- tidy=TRUE----------------------------------------------------------
names(GGE1)

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Biplot")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Selected Environment", selectedE = "OA93")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Selected Genotype", selectedG = "Kat")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Relationship Among Environments")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Comparison of Genotype", selectedG1 = "Kat", selectedG2 = "Cas")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Which Won Where/What")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Discrimination vs. representativeness")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Ranking Environments")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Ranking Genotypes")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
GGEPlot(GGE1, type="Mean vs. Stability")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan,genotype="gen",environment="env", response="yield", type = "AMMI")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan, genotype="gen",environment="env", response="yield", type = "rAMMI")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan, genotype="gen",environment="env", response="yield", type = "hAMMI")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan, genotype="gen",environment="env", response="yield", type = "gAMMI")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan, genotype="gen",environment="env", response="yield", type = "lAMMI")

## ----fig.height=5.5, fig.width=5.5, fig.align='center', tidy=TRUE--------
rAMMI(dat_yan, genotype="gen",environment="env", response="yield", type = "ppAMMI")

## ---- tidy=TRUE----------------------------------------------------------
# generates missing data
dat_yan[1,3]<-NA
dat_yan[3,3]<-NA
dat_yan[2,3]<-NA

## ---- tidy=TRUE----------------------------------------------------------
imputation(dat_yan, PC.nb=2, genotype="gen",environment="env", response="yield", type="EM-AMMI")

## ---- tidy=TRUE----------------------------------------------------------
imputation(dat_yan, PC.nb=1, genotype="gen",environment="env", response="yield", type="EM-AMMI")

## ---- tidy=TRUE----------------------------------------------------------
imputation(dat_yan, genotype="gen",environment="env", response="yield", type="EM-SVD")

## ---- tidy=TRUE----------------------------------------------------------
imputation(dat_yan, genotype="gen",environment="env", response="yield", type="WGabriel")

## ---- tidy=TRUE----------------------------------------------------------
imputation(dat_yan, genotype="gen",environment="env", response="yield", type="EM-PCA")

## ---- echo = FALSE-------------------------------------------------------
sessionInfo()

