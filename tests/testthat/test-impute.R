# library(testthat)


# Tests for impute function

library(geneticae)
# Data without replication
library(agridat)
data(yan.winterwheat)

dat2 <- yan.winterwheat
dat2[1,3]<-NA
dat2[3,3]<-NA
dat2[2,3]<-NA


test_that("Several tests for impute function", {

  expect_error(imputation(yan.winterwheat, genotype="gen",environment="env", response="yield", type="EM-AMMI"),
               "There are not missing data in input data frame")

  expect_error(imputation( genotype="gen",environment="env", response="yield", type="EM-AMMI"),
               "Need to provide Data data frame")


  expect_equal(imputation(dat2, genotype="gen",environment="env", response="yield", type="EM-AMMI"),
               imputation(dat2))


  # data types correct
  expect_that(dat2, is_a('data.frame'))

})
