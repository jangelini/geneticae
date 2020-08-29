
# Tests for GGE_Model function
library(geneticae)
# Data without replication
data(yan.winterwheat)
dat <- yan.winterwheat
GGE1 <- GGEmodel(dat, genotype="gen",environment="env", response="yield", rep=NULL, centering = "tester")

dat2 <- dat
dat2[1,3]<-NA
dat2[3,3]<-NA
dat2[2,3]<-NA


test_that("Several tests for GGE_Model function", {

  expect_that(GGEmodel(dat,genotype="gen",environment="env", response="yield", centering = "tester"), is_a("GGEModel") )

  expect_error(GGEmodel(genotype="gen",environment="env", response="yield", centering = "tester"),
               "Need to provide Data data frame")

  expect_equal(GGEmodel(dat,genotype="gen",environment="env", response="yield", centering = "tester"),
              GGEmodel(dat))

  expect_identical(dat, na.omit(dat))

  expect_error(GGEmodel(dat2,genotype="gen",environment="env", response="yield", centering = "tester"),
               "Missing data in input data frame, run the imputation function first to complete the data set")

  # data types correct
  expect_that(dat, is_a('data.frame'))

})




