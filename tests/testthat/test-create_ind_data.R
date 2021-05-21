library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("creates right data",{
  set.seed(1234)
  test_data<-create_ind_data(N=10000, beta2=2, beta1=1)

  expect_equal(nrow(test_data),
               10000,
               label = "right number of rows")

  expect_type(test_data,
              "list")

  expect_snapshot_output(test_data)
})
