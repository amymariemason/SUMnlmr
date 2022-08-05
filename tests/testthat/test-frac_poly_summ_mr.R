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
