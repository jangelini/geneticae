
# Tests for GGE_Plot function
# Data without replication
library(geneticae)
# Data without replication
data(yan.winterwheat)
dat <- yan.winterwheat
GGE1 <- GGEmodel(dat, genotype="gen",environment="env", response="yield", rep=NULL, centering = "tester")

GGE2 <- GGEmodel(dat, genotype="gen",environment="env", response="yield", rep=NULL, centering = "none")


dat2 <- dat
dat2[1,3]<-NA
dat2[3,3]<-NA
dat2[2,3]<-NA


test_that("multiplication works", {
  expect_that(GGE1, is_a("GGEModel") )

  # expect_error(GGEPlot(GGE2), "GGEPlot is not compatible with GGE models produced without centering" )


})
