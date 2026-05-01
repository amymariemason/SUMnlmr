library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("processes summarised data correctly",{
  set.seed(1234)
  test_data<-LDL_CAD
  model <-with(LDL_CAD, piecewise_summ_mr(by=by, bx=bx, byse=byse, bxse=bxse,
                                          xmean=xmean, xmin=xmin, xmax=xmax,
                                          ci="bootstrap_se",
                                          nboot=1000,
                                          fig=TRUE,
                                          family="gaussian",
                                          ci_fig="ribbon",
                                          average.exposure.associations = TRUE)
  )

  expect_snapshot_output(summary(model))
  expect_snapshot_output(model$figure)

  model2 <-with(LDL_CAD, piecewise_summ_mr(by=by, bx=bx, byse=byse, bxse=bxse,
                                           xmean=xmean, xmin=xmin, xmax=xmax,
                                           ci="model_se",
                                           fig=TRUE,
                                           family="binomial",
                                           ci_fig="point",
                                           average.exposure.associations = TRUE)
  )

  expect_snapshot_output(summary(model2))
  expect_snapshot_output(model2$figure)

})

test_that("accepts create_nlmr_summary-style input directly", {
  set.seed(1234)
  test_data<-generated_data
  summ_data <- create_nlmr_summary(
    y = test_data$y.bin,
    x = test_data$X,
    g = test_data$g,
    covar = NULL,
    family = "binomial",
    q = 10,
    strata_method = "residual"
  )

  model <- piecewise_summ_mr(
    summ= summ_data,
    family = "binomial",
  )

  expect_s3_class(model, "piecewise_summ_mr")
})
