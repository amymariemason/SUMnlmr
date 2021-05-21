library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("processes summarised data correctly",{
  set.seed(1234)
  test_data<-LDL_CAD
  model <-with(LDL_CAD, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax,
                                          ci="bootstrap_se",
                                          nboot=1000,
                                          fig=TRUE,
                                          family="gaussian",
                                          ci_fig="ribbon")
  )

  expect_snapshot_output(summary(model))
  expect_snapshot_output(model$figure)

  model2 <-with(LDL_CAD, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax,
                                           ci="model_se",
                                           fig=TRUE,
                                           family="binomial",
                                           ci_fig="point")
  )

  expect_snapshot_output(summary(model2))
  expect_snapshot_output(model2$figure)

})
