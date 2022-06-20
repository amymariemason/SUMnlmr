library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("creates correct summary",{
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)
  expect_snapshot_warning(getGRvalues(X = round(test_data$X),
                             Zstratum=floor((rank(test_data$g,
                                                  ties.method = "random")-1)/10)+1))

  expect_snapshot_output(getGRvalues(X = round(test_data$X),
                                     Zstratum=floor((rank(test_data$g,
                                          ties.method = "random")-1)/4)+1))
})

#first should produce GR warning, second should be fine
