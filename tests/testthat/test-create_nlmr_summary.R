library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("creates correct summary",{
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)
  expect_snapshot_output(create_nlmr_summary(y = test_data$log.Y,
                                                        x = test_data$X,
                                                        g = test_data$g,
                                                        covar = NULL,
                                                        family = "gaussian",
                                                        q = 10,
                                                        strata_method="residual")$summary)
})

test_that("creates correct summary2",{
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)

  expect_snapshot_output(create_nlmr_summary(y = test_data$y.bin,
                                             x = test_data$X,
                                             g = test_data$g,
                                             covar = as.matrix(test_data$linear.Y, ncol=1),
                                             family = "binomial",
                                             controlsonly = TRUE,
                                             q = 10,
                                             strata_method="residual")$summary)

  expect_snapshot_output(create_nlmr_summary(y = test_data$log.Y,
                                             x = test_data$X,
                                             g = test_data$g,
                                             covar = NULL,
                                             family = "gaussian",
                                             q = 10,
                                             strata_method = "ranked", seed=NA
                                            )$summary)

  expect_snapshot_output(create_nlmr_summary(y = test_data$log.Y,
                                             x = round(test_data$X),
                                             g = test_data$g,
                                             covar = NULL,
                                             family = "gaussian",
                                             q = 4,
                                             strata_method = "ranked",
                                             extra_statistics = TRUE,
                                             seed=NA)$strata_statistics)

  expect_snapshot_output(create_nlmr_summary(y = test_data$log.Y,
                                             x = round(test_data$X),
                                             g = test_data$g,
                                             covar = NULL,
                                             family = "gaussian",
                                             q = 4,
                                             strata_method = "ranked",
                                             report_GR = TRUE,
                                             seed=NA
                                            )$GR_results)
})

test_that("throws errors", {
  set.seed(1234)
  test_data<-generated_data
 expect_error(create_nlmr_summary(y = test_data$log.Y,
                                  x = test_data$X,
                                  g = test_data$g,
                                  covar = NULL,
                                  family = "smocks",
                                  q = 10,
                                  strata_method="residual"))



  expect_error(create_nlmr_summary(y = 1,
                                   x = test_data$X,
                                   g = test_data$g,
                                   covar = NULL,
                                   family = "gaussian",
                                   q = 10,
                                   strata_method="residual"))


  expect_snapshot_warning(create_nlmr_summary(y = test_data$log.Y,
                                             x = round(test_data$X),
                                             g = test_data$g,
                                             covar = NULL,
                                             family = "gaussian",
                                             q = 20,
                                             strata_method = "ranked"))

})
