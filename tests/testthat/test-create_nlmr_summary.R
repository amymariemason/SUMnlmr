library(testthat)
local_edition(3)
library(SUMnlmr)


test_that("creates correct summary",{
  set.seed(1234)
  test_data<-generated_data
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


test_that("supports custom strata as input", {
  set.seed(1234)
  test_data <- generated_data
  e_cov <- as.matrix(test_data$u)

  summ_strata <- create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    x_strata = rep(1:20, length.out = 10000),
    seed = 1234
  )

  expect_equal(nrow(summ_strata$summary), 20)

  expect_snapshot_output(create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    x_strata = rep(1:20, length.out = 10000),
    seed = 1234
  ))


}
)

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
    x_corrected = test_data$u,
    seed = 1234
  )

  expect_length(summ_xs$summary, 7)

  expect_snapshot_output(create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    x_corrected = test_data$X - test_data$u,
    seed = 1234
  ))


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
    gxe_covar = e_cov,
    gxe_interaction = e_interaction,
    extra_statistics = TRUE,
    report_het =TRUE,
    seed = 1234
  )

  expect_true("strata_statistics" %in% names(summ_method2))
  expect_true("Heterogeneity_results" %in% names(summ_method2))
  expect_equal(nrow(summ_method2$summary), 10)

  expect_snapshot_output(create_nlmr_summary(
    y = test_data$log.Y,
    x = test_data$X,
    g = test_data$g,
    q = 10,
    family = "gaussian",
    strata_method = "ranked",
    gxe_covar = e_cov,
    gxe_interaction = e_interaction,
    seed = 1234
  ))

})


test_that("gxe inputs run for gaussian, binomial and coxph", {
  set.seed(1234)
  test_data <- generated_data
  test_data$time <- test_data$log.Y + 10
  y_surv <- survival::Surv(time = test_data$time, event = test_data$y.bin)

  e_cov <- as.matrix(cbind(test_data$u, test_data$linear.Y))
  covar_mat <- as.matrix(cbind(test_data$linear.Y, test_data$quadratic.Y))
  gxe_term <- as.matrix(test_data$g * test_data$u)
  gxe_covar <- as.matrix(test_data$linear.Y)

  expect_snapshot_output(
    create_nlmr_summary(
      y = test_data$log.Y,
      x = test_data$X,
      g = test_data$g,
      covar = covar_mat,
      family = "gaussian",
      q = 6,
      strata_method = "ranked",
      gxe_interaction = gxe_term,
      gxe_covar= gxe_covar,
      seed = 1234
    )$summary
  )

  expect_snapshot_output(
    create_nlmr_summary(
      y = test_data$y.bin,
      x = test_data$X,
      g = test_data$g,
      covar = covar_mat,
      family = "binomial",
      controlsonly = TRUE,
      q = 6,
      strata_method = "ranked",
      gxe_covar = gxe_covar,
      gxe_interaction = e_cov,
      seed = 1234
    )$summary
  )

  expect_snapshot_output(
    create_nlmr_summary(
      y = y_surv,
      x = test_data$X,
      g = test_data$g,
      covar = covar_mat,
      family = "coxph",
      q = 6,
      strata_method = "ranked",
      gxe_covar = e_cov,
      gxe_interaction = gxe_covar,
      seed = 1234
    )$summary
  )
})

