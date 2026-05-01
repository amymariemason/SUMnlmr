library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("processes summarised data correctly",{

  test_data<-LDL_CAD
  model<- with(LDL_CAD, frac_poly_summ_mr(bx=bx,
                                                    by=by,
                                                    bxse=bxse,
                                                    byse=byse,
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                          average.exposure.associations = TRUE)
  )

  expect_snapshot_output(summary(model))
  expect_snapshot_output(model$figure)

  model2<- with(LDL_CAD, frac_poly_summ_mr(bx=bx,
                                          by=by,
                                          bxse=bxse,
                                          byse=byse,
                                          xmean=xmean,
                                          family="gaussian",
                                          fig=TRUE,
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

  model <- frac_poly_summ_mr(
    summ= summ_data,
    family = "binomial",
  )

  expect_s3_class(model, "frac_poly_mr")
})
