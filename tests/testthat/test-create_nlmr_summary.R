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

test_that("supports custom stratification inputs", {
  set.seed(1234)
  test_data <- generated_data
  e_cov <- as.matrix(test_data$u)

  summ_xs <- create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    xs = test_data$X + test_data$u,
    seed = 1234
  )

  expect_length(summ_xs$strata, length(test_data$X))
  expect_true(all(!is.na(summ_xs$strata)))

  expect_snapshot_output(create_nlmr_summary(create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    xs = test_data$X + test_data$u,
    seed = 1234
  )))


}
)


test_that("creates custom stratification inputs", {

  set.seed(1234)
  test_data <- generated_data
  e_cov <- as.matrix(test_data$u)
  e_interaction <- as.matrix(test_data$u^2)

  summ_method2 <- create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    xs_covar = e_cov,
    xs_interaction = e_interaction,
    include_standard_output = TRUE,
    seed = 1234
  )

  expect_true("summary_standard" %in% names(summ_method2))
  expect_true("strata_standard" %in% names(summ_method2))
  expect_equal(nrow(summ_method2$summary_standard), 10)
  expect_length(summ_method2$strata_standard, length(test_data$X))


  expect_snapshot_output(create_nlmr_summary(create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    xs_covar = e_cov,
    xs_interaction = e_interaction,
    include_standard_output = TRUE,
    seed = 1234
  )))

})

