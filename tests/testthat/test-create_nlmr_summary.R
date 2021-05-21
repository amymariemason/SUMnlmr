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
                                                        q = 10)$summary)

  expect_snapshot_output(create_nlmr_summary(y = test_data$y.bin,
                                             x = test_data$X,
                                             g = test_data$g,
                                             covar = as.matrix(test_data$linear.Y, ncol=1),
                                             family = "binomial",
                                             q = 10)$summary)
})

test_that("throws errors", {

 expect_error(create_nlmr_summary(y = test_data$log.Y,
                                  x = test_data$X,
                                  g = test_data$g,
                                  covar = NULL,
                                  family = "smocks",
                                  q = 10))



  expect_error(create_nlmr_summary(y = 1,
                                   x = test_data$X,
                                   g = test_data$g,
                                   covar = NULL,
                                   family = "gaussian",
                                   q = 10))


  expect_error(create_nlmr_summary(y = test_data$log.Y,
                                   x = test_data$X,
                                   g = test_data$g,
                                   covar = NULL,
                                   family = "gaussian",
                                   q = 1))

})
