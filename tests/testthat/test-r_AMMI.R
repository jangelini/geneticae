
# Tests for r_AMMI function

library(geneticae)
# Data without replication
data(yan.winterwheat)

dat2 <- yan.winterwheat
dat2[1,3]<-NA
dat2[3,3]<-NA
dat2[2,3]<-NA




test_that("Several tests for r_AMMI function", {

  expect_error(rAMMI(dat2, genotype="gen", environment="env", response="Y", rep=NULL),
               "Missing data in input data frame, run the imputation function first to complete the data set")

  expect_error(rAMMI(genotype="gen", environment="env", response="Y", rep=NULL),
               "Need to provide Data data frame")

  # data types correct
  expect_that(yan.winterwheat, is_a('data.frame'))

})
