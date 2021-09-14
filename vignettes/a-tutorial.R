## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(warning = F, message = F, out.width = "60%")

## ---- eval=F------------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("jangelini/geneticae")

## -----------------------------------------------------------------------------
library(geneticae)

## -----------------------------------------------------------------------------
data(yan.winterwheat)
head(yan.winterwheat)

## -----------------------------------------------------------------------------
data(plrv)
head(plrv)

## ---- fig.align='center', fig.cap='Figure 1: GE biplot.'----------------------
rAMMI(yan.winterwheat, genotype = "gen", environment = "env", 
      response = "yield", type = "AMMI", footnote = F, titles = F)

## -----------------------------------------------------------------------------
GGE1 <- GGEmodel(yan.winterwheat, genotype = "gen", environment = "env", 
                 response = "yield", centering = "tester")

## ---- fig.align='center', fig.cap='Figure 2: basic biplot. 78% of the total variability of G + GE is explained.'----
GGEPlot(GGE1, type="Biplot", footnote = F, titles = F)

## ---- fig.align='center', fig.cap='Figure 3: comparison of cultivare performance in a selected environment (OA93).'----
GGEPlot(GGE1, type="Selected Environment", selectedE = "OA93", 
        footnote = F, titles = F)

## ----  fig.align='center', fig.cap='Figure 4: comparison of the performance of cultivar Luc in different environments.'----
GGEPlot(GGE1, type="Selected Genotype", selectedG = "Luc", 
        footnote = F, titles = F)

## ----  fig.align='center', fig.cap='Figure 5: comparison of the cultivars Kat and Cas.'----
GGEPlot(GGE1, type = "Comparison of Genotype", 
        selectedG1 = "Kat", selectedG2 = "Cas", 
        footnote = F, titles = F, axis_expand = 1.5)

## ----  fig.align='center', fig.cap='Figure 6: polygon view of the GGE biplot, showing which cultivars presented highest yiel in each environment.'----
GGEPlot(GGE1, type="Which Won Where/What", footnote = F, titles = F, axis_expand = 1.5)


## ----  fig.align='center', fig.cap='Figure 7: average environment view of the GGE biplot based on genotype-focused scaling, showing mean yield and stability of genotypes.'----
data <- yan.winterwheat[yan.winterwheat$env %in% c("BH93", "EA93","HW93", "ID93",
                                                   "NN93", "RN93", "WP93"), ]
GGE2 <- GGEmodel(data, genotype = "gen", environment = "env", 
                 response = "yield", SVP = "row")
GGEPlot(GGE2, type = "Mean vs. Stability", footnote = F, titles = F)

## ----  fig.align='center', fig.cap='Figure 8: Classification of genotypes with respect to the ideal genotype.', warning=FALSE----
GGEPlot(GGE1, type = "Ranking Genotypes", footnote = F, titles = F)

## ----  fig.align='center', fig.cap='Figure 9: Relationship between environments.'----
GGE3 <- GGEmodel(data, genotype = "gen", environment = "env", 
                 response = "yield", SVP = "column")
GGEPlot(GGE3, type = "Relationship Among Environments", footnote = F, titles = F)

## ----  fig.align='center', fig.cap='Figure 10: classification of environments with respect to the ideal environment.'----
GGEPlot(GGE1, type = "Ranking Environments", footnote = F, titles = F)

## ----  fig.align='center'-----------------------------------------------------
GGEPlot(GGE1, type="Discrimination vs. representativeness", footnote = F, titles = F)

## -----------------------------------------------------------------------------
# Generating missing data
yan.winterwheat[1,3] <- NA
yan.winterwheat[3,3] <- NA
yan.winterwheat[2,3] <- NA

## -----------------------------------------------------------------------------
imputation(yan.winterwheat, nPC = 2, genotype = "gen", environment = "env", 
           response = "yield", type = "EM-AMMI")

## ---- echo = FALSE------------------------------------------------------------
sessionInfo()

