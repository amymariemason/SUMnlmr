library(testthat)
local_edition(3)
library(SUMnlmr)
library(survival)

test_that("throws errors if not Surv format", {
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)
  test_data$time<-test_data$log.Y+10
  new.y<- Surv(time = test_data$time, event = test_data$y.bin)
  expect_error(create_nlmr_summary(y = y.bin,
                                   x = test_data$X,
                                   g = test_data$g,
                                   covar = as.matrix(test_data$linear.Y, ncol=1),
                                   family = "coxph",
                                   q = 10,
                                   strata_method="residual"))
})

test_that("throws errors if coxph & controls only", {
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)
  test_data$time<-test_data$log.Y+10
  new.y<- Surv(time = test_data$time, event = test_data$y.bin)
  expect_error(create_nlmr_summary(y = new.y,
                                   x = test_data$X,
                                   g = test_data$g,
                                   covar = as.matrix(test_data$linear.Y, ncol=1),
                                   family = "coxph",
                                   controlsonly = TRUE,
                                   q = 10,
                                   strata_method="residual"))
})


test_that("coxph method running", {
  set.seed(1234)
  test_data<-generated_data
  test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)
  test_data$time<-test_data$log.Y+10
  new.y<- Surv(time = test_data$time, event = test_data$y.bin)

  expect_snapshot_output(create_nlmr_summary(y = new.y,
                                             x = test_data$X,
                                             g = test_data$g,
                                             covar = as.matrix(test_data$linear.Y, ncol=1),
                                             family = "coxph",
                                             q = 10,
                                             strata_method="residual")$summary)

})
