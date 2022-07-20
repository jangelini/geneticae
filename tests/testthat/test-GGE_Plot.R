
# Tests for GGE_Plot function
# Data without replication
library(geneticae)
# Data without replication
library(agridat)
data(yan.winterwheat)

GGE1 <- GGEmodel(yan.winterwheat, genotype="gen",environment="env", response="yield", rep=NULL)

GGE2 <- GGEmodel(yan.winterwheat, genotype="gen",environment="env", response="yield", rep=NULL)



test_that("multiplication works", {
  expect_that(GGE1, is_a("GGEModel") )

  expect_error(GGEPlot(GGE1,type="Selected Genotype",selectedG="none"), "The genotype selected is not in list of genotype labels" )
  expect_error(GGEPlot(GGE1,type="Comparison of Genotype", selectedG1="none", selectedG2="Fun"), "none The genotype selected is not in list of genotype labels" )
  expect_error(GGEPlot(GGE1,type="Comparison of Genotype", selectedG1="none", selectedG2="none"), "none The genotype selected is not in list of genotype labels")
  expect_error(GGEPlot(GGE1,type="Comparison of Genotype", selectedG1="Fun", selectedG2="Fun"), "Cannot compare the same genotype to itself")
  expect_error(GGEPlot(GGE1,type="Selected Environment", selectedE="none"), "The environment selected is not in list of environment labels")
})
