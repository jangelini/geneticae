# library(testthat)


# Tests for impute function

library(geneticae)
# Data without replication
data(yan.winterwheat)
dat <- yan.winterwheat

dat2 <- dat
dat2[1,3]<-NA
dat2[3,3]<-NA
dat2[2,3]<-NA


test_that("Several tests for impute function", {

  expect_error(imputation(dat, genotype="gen",environment="env", response="yield", type="EM-SVD"),
               "There are not missing data in input data frame")

  expect_error(imputation( genotype="gen",environment="env", response="yield", type="EM-SVD"),
               "Need to provide Data data frame")


  expect_equal(imputation(dat2, genotype="gen",environment="env", response="yield", type="EM-AMMI"),
               imputation(dat2))


  # data types correct
  expect_that(dat2, is_a('data.frame'))

})
