context("Matrix test")

## Check that conditions are a vector
  mat <- matrix(NA, nrow=100, ncol=10)
  conds <- matrix(NA, nrow=10, ncol=1)
  expect_error(HTSFilter(mat, conds))

test_that("data.frame test", {
  expect_that(10, equals(10))
  expect_that(11, equals(11))
})
