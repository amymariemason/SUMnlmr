library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("creates correct warning 1",{
  set.seed(1234)
  test_data<-generated_data
  expect_snapshot_warning(getGRvalues(X = round(test_data$X),
                             Zstratum=floor((rank(test_data$g,
                                                  ties.method = "random")-1)/20)+1))
})

test_that("creates correct summary 2",{
  set.seed(1234)
  test_data<-generated_data
  expect_snapshot_output(getGRvalues(X = round(test_data$X),
                                     Zstratum=floor((rank(test_data$g,
                                          ties.method = "random")-1)/4)+1))
})

#first should produce GR warning, second should be fine
