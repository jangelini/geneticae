# Tests for geneticae Package

data(Ontario)

results <- GGEmodel(geneticae)

test_that("Several tests", {
    expect_error(GGEmodel(Ontario, rep=2),
                 'rep is logical', fixed = TRUE)
})
