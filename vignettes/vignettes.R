## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("jangelini/geneticae")

## ------------------------------------------------------------------------
library(geneticae)
library(reshape2)
data(yan.winterwheat)
dat1 <- yan.winterwheat
dat <-t(round(acast(dat1, env~gen, value.var='yield'),2))

GGE1<-GGEmodel(dat, centering="tester", rep=FALSE)
names(GGE1)
GGE1$varexpl
GGE1$centering
GGE1$scaling
GGE1$SVP

## ----fig.height=5.5, fig.width=5.5---------------------------------------
GGEPlot(GGE1, type="Which Won Where/What")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
GGEPlot(GGE1, type="Mean vs. Stability")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
rAMMI(dat1,type='AMMI', rep=FALSE)

## ----fig.height=5.5, fig.width=5.5---------------------------------------
rAMMI(dat1,type='rAMMI', rep=FALSE)

## ----fig.height=5.5, fig.width=5.5---------------------------------------
# generates missing data
dat[1,2]<-NA
dat[3,4]<-NA
dat[2,7]<-NA

## ----fig.height=5.5, fig.width=5.5---------------------------------------
imputation(dat, rep=FALSE, type="EM-SVD")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
imputation(dat, rep=FALSE, type="EM-AMMI")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
imputation(dat, rep=FALSE, type="Gabriel")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
imputation(dat, rep=FALSE, type="WGabriel")

## ----fig.height=5.5, fig.width=5.5---------------------------------------
imputation(dat, rep=FALSE, type="EM-PCA")

## ---- echo = FALSE-------------------------------------------------------
sessionInfo()

