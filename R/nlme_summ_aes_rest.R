if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "y", "yest"))
}
#' Instrumental variable analysis using fractional polynomials based on summary
#' data
#'
#' @description frac_poly_summ_mr performs instumental variable analysis by
#' fitting fractional polynomial models to localised average causal effects
#' using meta-regression.
#'
#' Please note that if you provide synthetic data that are too regular (eg the
#' by associations are all zero or the xmean values are exactly 1, 2, 3, ...),
#' the function may error as several fractional polynomials provide the same
#' fit.
#'
#' @param by vector of gene-outcome associations.
#' @param bx vector of gene-exposure associations.
#' @param byse vector of standard errors of gene-outcome associations.
#' @param bxse vector of standard errors of gene-exposure associations.
#' @param xmean average value of the exposure in each stratum (or whatever
#' summary of the exposure level in the stratum is desired).
#' @param method meta-regression method parsed to the rma package. The default
#' is fixed-effects ('FE').
#' @param d fractional polynomial degree. The default is degree 1. The other
#' options are: 1, 2, or 'both'.
#' @param powers fractional polynomial powers to test.
#' @param pd p-value cut-off for choosing the best-fitting fractional polynomial
#'  of degree 2 over the best-fitting fractional polynomial degree 1. This
#'  option is only used if d='both'. The default is 0.05.
#' @param ci the type of 95\\% confidence interval. There are three options:
#' (i) using the model standard errors ('model_se'), (ii) using bootstrap
#' standard errors ('bootstrap_se'), (iii) using bootstrap percentile
#' confidence intervals ('bootstrap_per'). The default is the model standard
#' errors.
#' @param nboot the number of bootstrap replications (if required). The default
#'  is 100 replications.
#' @param fig a logical statement as to whether the user wants the results
#' displayed in a figure. The default is false.
#' @param family a character string named either \code{'gaussian'} (for
#' continuous data) or binomial (for binary data) family function. This only
#' affects the plotting function - whether the y-axis is log-transformed of not.
#' @param offset offset on the x-axis (default is zero).
#' @param pref_x the prefix/label for the x-axis. The default is 'x'.
#' @param pref_y the prefix/label for the y-axis. The default is 'y'.
#' @param ref the reference point for the figure. This is the value of the
#' function that represents the expected difference in the outcome compared
#' with this reference value when the exposure is set to different values.
#' If \code{ref = NA} (the default option), then it is set to the mean of x.
#' @param ci_type the type of confidence interval to be displayed on the graph.
#' The default is 'overall' where confidence intervals are presented as bands
#' across the range of x. The alternative option is 'quantile' where the
#' confidence intervals are presented as error bars at the mean in each
#' quantile of x.
#' @param breaks breaks on the y-axis of the figure.
#' @param xlim_lower lower limit for the x-axis of the figure.
#' @param xlim_upper upper limit for the x-axis of the figure.
#' @param ylim_lower lower limit for the y-axis of the figure.
#' @param ylim_upper upper limit for the y-axis of the figure.
#' @return model the model specifications. The first column is the number of
#' quantiles (q); the second column is the position used to relate x to the LACE
#'  in each quantiles (xpos); the third column is the type of confidence
#'  interval constructed (ci); the fourth column is the number of bootstrap
#'  replications performed (nboot).
#' @return powers the powers of the chosen polynomial.
#' @return coefficients the regression estimates. The first column is the
#' regression coefficients (beta); the second column is the standard errors of
#' regression coefficients (se); the third column is the lower confidence
#' interval (lci); the fourth column is the upper confidence interval (uci);
#' the fifth column is the p-value (pval).
#' @return lace the localised average causal effect estimate in each quantile.
#' The first column is the regression coefficients (beta); the second column is
#' the standard errors of regression coefficients (se); the third column is the
#'  lower confidence interval (lci); the fourth column is the upper confidence
#'   interval (uci); the fifth column is the p-value (pval).
#' @return xcoef the association between the instrument and the exposure in each
#'  quantile. The first column is the regression coefficients (beta); the second
#'   column is the standard errors of regression coefficients (se).
#' @return p_tests the p-value of the non-linearity tests. The first column is
#' the p-value of the test between the fractional polynomial degrees (fp_d1_d2);
#' the second column is the p-value from the fractional polynomial non-linearity
#'  test (fp); the third column is the p-value from the quadratic test (quad);
#'  the fourth column is the p-value from the Cochran Q test (Q).
#' @return p_heterogeneity the p-value of heterogeneity.
#' The first column is the p-value of the Cochran Q heterogeneity test (Q);
#' the second column is the p-value from the trend test (trend).
#' @return figure ggplot command to produce a figure.
#' @author Stephen Burgess <sb452@medschl.cam.ac.uk>, leading heavily on
#' James R Staley <js16174@bristol.ac.uk>
#' @import ggplot2
#' @importFrom matrixStats rowQuantiles
#' @importFrom metafor rma
#' @importFrom metafor rma.uni
#' @export

## Matt Arnold <mga37@medschl.cam.ac.uk> changes - add the powers option -
## fix in line 204 to correct a typo; add $x to give
## log(plot_data$x) in the second from last term

frac_poly_summ_mr <- function(by, bx, byse, bxse, xmean, method = "FE", d = 1,
                              powers = c(0, -3, -2, -1.5, -1, -0.5, 1, 2),
                              pd = 0.05, ci = "model_se", nboot = 100,
                              fig = FALSE, family = "binomial", offset = 0,
                              pref_x = "x", pref_y = "y", ref = NA,
                              ci_type = "overall", breaks = NULL,
                              ylim_lower = NA, ylim_upper = NA,
                              xlim_lower = NA, xlim_upper = NA) {
    ##### Error messages #####
    if (!(d == 1 | d == 2 | d == "both"))
        stop("the degree has to equal 1, 2 or \"both\"")
    if (!(ci == "model_se" | ci == "bootstrap_se" | ci == "bootstrap_per")) {
        stop("the confidence intervals must be one of \"model_se\",
         \"bootstrap_se\" and \"bootstrap_per\"")
    }
    #### Start ###
    frac_coef <- by
    frac_se <- byse
    xcoef_sub <- bx
    xcoef_sub_se <- bxse
    xcoef <- sum(bx * (bxse^-2))/sum(bxse^-2)
    q <- length(by)

    ##### Test of IV-exposure assumption #####
    p_het <- 1 - stats::pchisq(rma(xcoef_sub,
                                            vi = (xcoef_sub_se)^2)$QE,
                                            df = (q - 1))
    p_het_trend <- rma.uni(xcoef_sub ~ xmean,
                                    vi = xcoef_sub_se^2,
                                    method = method)$pval[2]

    ##### Best-fitting fractional polynomial of degree 1 #####
    p <- NULL
    ML <- NULL
    j <- 1
    for (p1 in powers) {
        if (p1 == -1) {
            x1 <- xmean^p1
        } else {
            x1 <- (p1 + 1) * xmean^p1
        }
        p[j] <- p1
        ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1,
                                      vi = (frac_se/xcoef)^2,
                                      method = method))$fit.stats[1, 1]
        j <- j + 1
    }
    fits <- data.frame(p, ML)
    fits$max <- 0
    fits$max[fits$ML == max(fits$ML)] <- 1
    p_ML <- fits$p[fits$max == 1]

    ##### Best-fitting fractional polynomial of degree 2 #####
    if (d == 1 | d == 2 | d == "both") {
        powers1 <- powers
        powers2 <- powers
        p1 <- NULL
        p2 <- NULL
        ML <- NULL
        j <- 1
        for (p11 in powers1) {
            if (p11 == -1) {
                x1 <- xmean^p11
            } else {
                x1 <- (p11 + 1) * xmean^p11
            }
            for (p21 in powers2) {
                if (p11 == p21) {
                  if (p21 == -1) {
                    x2 <- 2 * (xmean^p21) * log(xmean)
                  } else {
                    x2 <- ((p21 + 1) * (xmean^p21) * log(xmean) + xmean^p21)
                  }
                } else {
                  if (p21 == -1) {
                    x2 <- xmean^p21
                  } else {
                    x2 <- (p21 + 1) * xmean^p21
                  }
                }
                p1[j] <- p11
                p2[j] <- p21
                cc <- try(rma(frac_coef/xcoef ~ -1 + x1 + x2,
                                       vi = (frac_se/xcoef)^2,
                                       method = method), silent = TRUE)
                if (methods::is(cc, "try-error") == T) {
                  ML[j] <- NA
                }
                if (methods::is(cc, "try-error") == F) {
                  ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1 + x2,
                                                vi = (frac_se/xcoef)^2,
                                                method = method))$fit.stats[1,
                    1]
                }
                j <- j + 1
            }
            powers2 <- powers2[-1]
        }
        fits <- data.frame(p1, p2, ML)
        fits$max <- 0
        fits$max[fits$ML == max(fits$ML, na.rm = T)] <- 1
        p1_ML <- fits$p1[fits$max == 1]
        p2_ML <- fits$p2[fits$max == 1]
    }

    ##### Best-fitting fractional polynomial of either degree 1 or degree 2 ###
    p_d1_d2 <- NA
    if (d == 1 | d == 2 | d == "both") {
        if (p_ML == -1) {
            x1 <- xmean^p_ML
        } else {
            x1 <- (p_ML + 1) * xmean^p_ML
        }
        best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1,
                                      vi = (frac_se/xcoef)^2,
                                      method = method)
        dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2, 1]
        if (p1_ML == -1) {
            x1 <- xmean^p1_ML
        } else {
            x1 <- (p1_ML + 1) * xmean^p1_ML
        }
        if (p1_ML == p2_ML) {
            if (p2_ML == -1) {
                x2 <- 2 * (xmean^p2_ML) * log(xmean)
            } else {
                x2 <- ((p2_ML + 1) * (xmean^p2_ML) * log(xmean) + xmean^p2_ML)
            }
        } else {
            if (p2_ML == -1) {
                x2 <- xmean^p2_ML
            } else {
                x2 <- (p2_ML + 1) * xmean^p2_ML
            }
        }
        best_fracp_d2 <- rma(frac_coef/xcoef ~ -1 + x1 + x2,
                                      vi = (frac_se/xcoef)^2,
                                      method = method)
        dev_best_fracp_d2 <- best_fracp_d2$fit.stats[2, 1]
        p_d1_d2 <- 1 - stats::pchisq((dev_best_fracp_d1 - dev_best_fracp_d2),
                                     df = 2)
        if (p_d1_d2 >= pd) {
            d1 <- 1
        } else {
            d1 <- 2
        }
        if (d == "both") {
            d <- d1
        }
    }

    ##### Model #####
    if (d == 1) {
        if (p_ML == -1) {
            x1 <- xmean^p_ML
        } else {
            x1 <- (p_ML + 1) * xmean^p_ML
        }
        model <- rma(frac_coef/xcoef ~ -1 + x1,
                              vi = (frac_se/xcoef)^2,
                              method = method)
    }
    if (d == 2) {
        if (p1_ML == -1) {
            x1 <- xmean^p1_ML
        } else {
            x1 <- (p1_ML + 1) * xmean^p1_ML
        }
        if (p1_ML == p2_ML) {
            if (p2_ML == -1) {
                x2 <- 2 * (xmean^p2_ML) * log(xmean)
            } else {
                x2 <- ((p2_ML + 1) * (xmean^p2_ML) * log(xmean) + xmean^p2_ML)
            }
        } else {
            if (p2_ML == -1) {
                x2 <- xmean^p2_ML
            } else {
                x2 <- (p2_ML + 1) * xmean^p2_ML
            }
        }
        model <- rma(frac_coef/xcoef ~ -1 + x1 + x2,
                              vi = (frac_se/xcoef)^2,
                              method = method)
    }

    ##### Bootstrap #####

    if (ci == "bootstrap_per" | ci == "bootstrap_se") {
        if (d == 1) {
            frac_coef_boot <- NULL
            for (i in 1:nboot) {
                frac_coef1 <- by + stats::rnorm(q, 0, byse)
                frac_se1 <- byse
                xmean1 <- xmean
                if (p_ML == -1) {
                  x111 <- xmean1^p_ML
                } else {
                  x111 <- (p_ML + 1) * xmean1^p_ML
                }
                mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111,
                               vi = (frac_se1/xcoef)^2, method = method)
                frac_coef_boot[i] <- mod$b[1]
            }
        }
        if (d == 2) {
            frac_coef_boot <- matrix(, nrow = nboot, ncol = 2)
            for (i in 1:nboot) {
                frac_coef1 <- by + stats::rnorm(q, 0, byse)
                frac_se1 <- byse
                xmean1 <- xmean
                if (p1_ML == -1) {
                  x111 <- xmean1^p1_ML
                } else {
                  x111 <- (p1_ML + 1) * xmean1^p1_ML
                }
                if (p1_ML == p2_ML) {
                  if (p2_ML == -1) {
                    x211 <- 2 * (xmean1^p2_ML) * log(xmean1)
                  } else {
                    x211 <- ((p2_ML + 1) * (xmean1^p2_ML) * log(xmean1)
                             + xmean1^p2_ML)
                  }
                } else {
                  if (p2_ML == -1) {
                    x211 <- xmean1^p2_ML
                  } else {
                    x211 <- (p2_ML + 1) * xmean1^p2_ML
                  }
                }
                mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111 + x211,
                               vi = (frac_se1/xcoef)^2, method = method)
                frac_coef_boot[i, 1] <- mod$b[1]
                frac_coef_boot[i, 2] <- mod$b[2]
            }
        }
    }

    ##### Fractional polynomial degree 1 test against linearity #####
    if (p_ML == -1) {
        x1 <- xmean^p_ML
    } else {
        x1 <- (p_ML + 1) * xmean^p_ML
    }
    linear <- rma(frac_coef/xcoef ~ 1, vi = (frac_se/xcoef)^2, method = method)
    dev_linear <- linear$fit.stats[2, 1]
    best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1, vi = (frac_se/xcoef)^2,
                         method = method)
    dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2, 1]
    p_fp <- 1 - stats::pchisq((dev_linear - dev_best_fracp_d1), df = 1)

    ##### Other tests #####
    p_quadratic <- rma(frac_coef/xcoef ~ xmean, vi = (frac_se/xcoef)^2,
                       method = method)$pval[2]
    p_q <- 1 - stats::pchisq(rma(frac_coef/xcoef, vi = (frac_se/xcoef)^2)$QE,
                             df = (q - 1))

    ##### Results #####
    beta <- as.numeric(model$b)
    if (ci == "model_se") {
        if (d == 1) {
            powers <- p_ML + 1
        }
        if (d == 2) {
            powers <- c(p1_ML, p2_ML)
            powers <- powers + 1
        }
        cov <- model$vb
        se <- model$se
        lci <- beta - 1.96 * se
        uci <- beta + 1.96 * se
        pval <- 2 * stats::pnorm(-abs(beta/se))
    }
    if (ci == "bootstrap_se") {
        if (d == 1) {
            powers <- p_ML + 1
            cov <- stats::var(frac_coef_boot)
            se <- sqrt(cov)
        }
        if (d == 2) {
            powers <- c(p1_ML, p2_ML)
            powers <- powers + 1
            cov <- cov(frac_coef_boot)
            se <- sqrt(diag(cov))
        }
        lci <- beta - 1.96 * se
        uci <- beta + 1.96 * se
        pval <- 2 * stats::pnorm(-abs(beta/se))
    }
    if (ci == "bootstrap_per") {
        if (d == 1) {
            powers <- p_ML + 1
            se <- NA
            lci <- stats::quantile(frac_coef_boot, probs = 0.025)
            uci <- stats::quantile(frac_coef_boot, probs = 0.975)
            pval <- NA
        }
        if (d == 2) {
            powers <- c(p1_ML, p2_ML)
            powers <- powers + 1
            se <- rep(NA, 2)
            lci <- NULL
            uci <- NULL
            pval <- NULL
            lci[1] <- stats::quantile(frac_coef_boot[, 1], probs = 0.025)
            lci[2] <- stats::quantile(frac_coef_boot[, 2], probs = 0.025)
            uci[1] <- stats::quantile(frac_coef_boot[, 1], probs = 0.975)
            uci[2] <- stats::quantile(frac_coef_boot[, 2], probs = 0.975)
            pval <- rep(NA, 2)
        }
    }
    lci <- as.numeric(lci)
    uci <- as.numeric(uci)
    if (ci == "model_se") {
        nboot <- NA
    }

    ##### Figure #####
    if (fig == TRUE) {
        if (ci_type == "overall") {
            if (is.na(ref)) {
                ref <- mean(xmean)
            }
            plot_data <- data.frame(x = stats::runif(10000, min(xmean),
                                                     max(xmean)))
            plot_data_1 <- data.frame(x = ref, y = 0)
            if (d == 1) {
                if (p_ML == -1) {
                  plot_data$yest <- beta * log(plot_data$x) - (beta * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) - log(ref))^2 *
                                            as.vector(cov))
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot) -
                      reprow(log(ref) %*% t(frac_coef_boot),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                } else {
                  plot_data$yest <- beta * plot_data$x^(p_ML + 1) -
                    beta * ref^(p_ML + 1)
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p_ML + 1) -
                                             ref^(p_ML + 1))^2 *  as.vector(cov))
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p_ML + 1) %*% t(frac_coef_boot) -
                      reprow(ref^(p_ML + 1) %*% t(frac_coef_boot),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
            }
            if (d == 2) {
                if (p1_ML == -1 & p2_ML == -1) {
                  plot_data$yest <- beta[1] * log(plot_data$x) +
                    beta[2] * log(plot_data$x) * log(plot_data$x) -
                    (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) -
                                          log(ref))^2 * cov[1, 1] +
                                          2 * (log(plot_data$x) -
                                          log(ref)) * (log(plot_data$x) *
                                                             log(plot_data$x) -
                                          log(ref) * log(ref)) * cov[1, 2] +
                                          (log(plot_data$x) * log(plot_data$x) -
                                             log(ref) * log(ref))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
                      log(plot_data$x) *
                      log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
                      reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                               log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML == -1 & p2_ML != -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * log(plot_data$x) +
                    beta[2] * plot_data$x^(p2_ML + 1) - (beta[1] * log(ref) +
                                                beta[2] * ref^(p2_ML +  1))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) -
                                          log(ref))^2 * cov[1, 1] +
                                          2 * (log(plot_data$x) -
                                          log(ref)) * (plot_data$x^(p2_ML + 1) -
                                          ref^(p2_ML + 1)) * cov[1, 2] +
                                          (plot_data$x^(p2_ML + 1) -
                                             ref^(p2_ML + 1))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
                      plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
                      reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                               ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML == -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * log(plot_data$x) -
                    (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                        ref^(p1_ML + 1))^2 * cov[1, 1] +
                                        2 * (plot_data$x^(p1_ML + 1) -
                                        ref^(p1_ML + 1)) * (log(plot_data$x) -
                                                        log(ref)) * cov[1, 2] +
                                        (log(plot_data$x) -
                                           log(ref))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                              log(ref) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML != -1 & p1_ML == p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                    (beta[1] * ref^(p1_ML + 1) +
                       beta[2] * ref^(p2_ML + 1) * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                   ref^(p1_ML + 1))^2 * cov[1, 1] +
                                  2 * (plot_data$x^(p1_ML + 1) -
                                  ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) *
                                      log(plot_data$x) -
                                      ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
                                  (plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                                     ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      (plot_data$x^(p2_ML + 1) *
                         log(plot_data$x)) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                        (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
                      n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML != -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * plot_data$x^(p2_ML + 1) -
                    (beta[1] * ref^(p1_ML + 1) +  beta[2] * ref^(p2_ML + 1))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                  ref^(p1_ML + 1))^2 * cov[1, 1] +
                                  2 * (plot_data$x^(p1_ML + 1) -
                                  ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) -
                                  ref^(p2_ML + 1)) * cov[1, 2] +
                                  (plot_data$x^(p2_ML + 1) -
                                  ref^(p2_ML + 1))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                        n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
            }
            plot_data$x <- plot_data$x + offset
            plot_data_1$x <- plot_data_1$x + offset
            ref <- ref + offset
            if (family != "binomial") {
                figure <- ggplot2::ggplot(plot_data, aes(x = x))
                figure <- figure +
                  geom_hline(aes(yintercept = 0), colour = "grey") +
                  geom_line(aes(y = yest), color = "black") +
                  geom_line(aes(y = lci),  color = "grey") +
                  geom_line(aes(y = uci), color = "grey") +
                  theme_bw() +
                  labs(x = pref_x, y = pref_y) +
                  theme(axis.title.x = element_text(vjust = 0.5, size = 20),
                        axis.title.y = element_text(vjust = 0.5, size = 20),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18)) +
                  geom_point(aes(x = x, y = y),
                             data = plot_data_1,
                             colour = "red",
                             size = 4) +
                  theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
                if (!is.null(breaks)) {
                  suppressMessages(figure <- figure +
                                     scale_y_continuous(breaks = breaks))
                }
            }
            if (family == "binomial") {
                plot_data$yest <- exp(plot_data$yest)
                plot_data$uci <- exp(plot_data$uci)
                plot_data$lci <- exp(plot_data$lci)
                plot_data_1$y <- exp(0)
                figure <- ggplot2::ggplot(plot_data, aes(x = x))
                figure <- figure +
                  geom_hline(aes(yintercept = 1), colour = "grey") +
                  geom_line(aes(y = yest), color = "black") +
                  geom_line(aes(y = lci),  color = "grey") +
                  geom_line(aes(y = uci), color = "grey") +
                  theme_bw() +
                  labs(x = pref_x, y = pref_y) +
                  theme(axis.title.x = element_text(vjust = 0.5, size = 20),
                        axis.title.y = element_text(vjust = 0.5, size = 20),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18)) +
                  geom_point(aes(x = x, y = y),
                             data = plot_data_1,
                             colour = "red",
                             size = 4) +
                  theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
                if (!is.null(breaks)) {
                  figure <- figure + scale_y_continuous(breaks = breaks)
                }
                figure <- figure + coord_trans(y = "log")
            }
        }
        if (ci_type == "quantile") {
            if (is.na(ref)) {
                ref <- mean(xmean)
            }
            xmean_ci <- xmean
            plot_data <- data.frame(x = c(ref, xmean_ci))
            plot_data_1 <- data.frame(x = stats::runif(10000, min(xmean),
                                                       max(xmean)))
            if (d == 1) {
                if (p_ML == -1) {
                  plot_data$yest <- beta * log(plot_data$x) -
                    (beta * log(ref))
                  plot_data_1$yest <- beta * log(plot_data_1$x) -
                    (beta * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) - log(ref))^2 *
                                            as.vector(cov))
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot) -
                      reprow(log(ref) %*% t(frac_coef_boot),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p_ML != -1) {
                  plot_data$yest <- beta * plot_data$x^(p_ML + 1) -
                    beta * ref^(p_ML + 1)
                  plot_data_1$yest <- beta * plot_data_1$x^(p_ML + 1) -
                    beta * ref^(p_ML + 1)
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p_ML + 1) -
                                             ref^(p_ML + 1))^2 *  as.vector(cov))
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p_ML + 1) %*% t(frac_coef_boot) -
                      reprow(ref^(p_ML + 1) %*% t(frac_coef_boot),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
            }
            if (d == 2) {
                if (p1_ML == -1 & p2_ML == -1) {
                  plot_data$yest <- beta[1] * log(plot_data$x) +
                    beta[2] * log(plot_data$x) * log(plot_data$x) -
                    (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
                  plot_data_1$yest <- beta[1] * log(plot_data_1$x) +
                    beta[2] * log(plot_data_1$x) * log(plot_data_1$x) -
                    (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) -
                              log(ref))^2 * cov[1, 1] +
                              2 * (log(plot_data$x) -
                              log(ref)) * (log(plot_data$x) * log(plot_data$x) -
                              log(ref) * log(ref)) * cov[1, 2] +
                              (log(plot_data$x) * log(plot_data$x) -
                              log(ref) * log(ref))^2 *  cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
                      log(plot_data$x) *
                        log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
                      reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                              log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML == -1 & p2_ML != -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * log(plot_data$x) +
                    beta[2] * plot_data$x^(p2_ML + 1) -
                    (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
                  plot_data_1$yest <- beta[1] * log(plot_data_1$x) +
                    beta[2] * plot_data_1$x^(p2_ML + 1) -
                    (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((log(plot_data$x) -
                                        log(ref))^2 * cov[1, 1] +
                                        2 * (log(plot_data$x) - log(ref)) *
                                          (plot_data$x^(p2_ML + 1) -
                                             ref^(p2_ML + 1)) * cov[1, 2] +
                                          (plot_data$x^(p2_ML + 1) -
                                             ref^(p2_ML + 1))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
                      plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
                      reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                               ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML == -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * log(plot_data$x) -
                    (beta[1] * ref^(p1_ML + 1) + beta[2] *  log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                             ref^(p1_ML + 1))^2 * cov[1, 1] +
                                            2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (log(plot_data) -
                                                        log(ref)) * cov[1, 2] +
                                                        (log(plot_data$x) -
                                                        log(ref))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                               log(ref) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML != -1 & p1_ML == p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                    (beta[1] * ref^(p1_ML + 1) +
                       beta[2] * ref^(p2_ML + 1) * log(ref))
                  plot_data_1$yest <- beta[1] * plot_data_1$x^(p1_ML + 1) +
                    beta[2] * plot_data_1$x^(p2_ML + 1) * log(plot_data_1$x) -
                    (beta[1] * ref^(p1_ML + 1) +
                       beta[2] * ref^(p2_ML + 1) * log(ref))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                ref^(p1_ML + 1))^2 * cov[1, 1] +
                                  2 * (plot_data$x^(p1_ML + 1) -
                                ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) *
                                                     log(plot_data$x) -
                                ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
                                  (plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                                     ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      (plot_data$x^(p2_ML + 1) *
                         log(plot_data$x)) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                        (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
                      n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
                if (p1_ML != -1 & p2_ML != -1 & p1_ML != p2_ML) {
                  plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
                    beta[2] * plot_data$x^(p2_ML + 1) -
                    (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
                  plot_data_1$yest <- beta[1] * plot_data_1$x^(p1_ML + 1) +
                    beta[2] * plot_data_1$x^(p2_ML + 1) -
                    (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
                  if (ci != "bootstrap_per") {
                    plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                             ref^(p1_ML + 1))^2 * cov[1, 1] +
                                            2 * (plot_data$x^(p1_ML + 1) -
                                  ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) -
                                                ref^(p2_ML + 1)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1))^2 * cov[2, 2])
                    plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
                    plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
                  } else {
                    boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                      plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
                      reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                               ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                             n = nrow(plot_data))
                    plot_data$lci <- rowQuantiles(boot, probs = 0.025)
                    plot_data$uci <- rowQuantiles(boot, probs = 0.975)
                  }
                }
            }
            highlight <- c("red", rep("black", (nrow(plot_data) - 1)))
            plot_data$x <- plot_data$x + offset
            plot_data_1$x <- plot_data_1$x + offset
            if (family != "binomial") {
                figure <- ggplot2::ggplot(plot_data, aes(x = x))
                figure <- figure +
                  geom_hline(aes(yintercept = 0), colour = "grey") +
                  geom_line(aes(x = x, y = yest), color = "black",
                            data = plot_data_1) +
                  geom_errorbar(mapping = aes(x = x, ymin = lci, ymax = uci),
                                color = "grey", width = 0.025) +
                  geom_point(aes(y = yest), color = highlight, size = 4) +
                  theme_bw() +
                  labs(x = pref_x, y = pref_y) +
                  theme(axis.title.x = element_text(vjust = 0.5, size = 20),
                        axis.title.y = element_text(vjust = 0.5, size = 20),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18)) +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
                if (!is.null(breaks)) {
                  figure <- figure + scale_y_continuous(breaks = breaks)
                }
            }
            if (family == "binomial") {
                plot_data$yest <- exp(plot_data$yest)
                plot_data$uci <- exp(plot_data$uci)
                plot_data$lci <- exp(plot_data$lci)
                figure <- ggplot2::ggplot(plot_data, aes(x = x))
                plot_data_1$yest <- exp(plot_data_1$yest)
                figure <- figure +
                  geom_hline(aes(yintercept = 1), colour = "grey") +
                  geom_line(aes(x = x, y = yest), color = "black",
                            data = plot_data_1) +
                  geom_errorbar(mapping = aes(x = x, ymin = lci, ymax = uci),
                                color = "grey", width = 0.025) +
                  geom_point(aes(y = yest), color = highlight, size = 4) +
                  theme_bw() +
                  labs(x = pref_x, y = pref_y) +
                  theme(axis.title.x = element_text(vjust = 0.5, size = 20),
                        axis.title.y = element_text(vjust = 0.5, size = 20),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18)) +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
                if (!is.null(breaks)) {
                  figure <- figure + scale_y_continuous(breaks = breaks)
                }
                figure <- figure + coord_trans(y = "log")
            }
        }
        figure <- figure + theme(panel.border = element_blank(),
                                 axis.line = element_line())
        if (!is.na(ylim_lower) | !is.na(ylim_upper)) {
            figure <- figure +
              scale_y_continuous(limits = c(ylim_lower, ylim_upper),
                                 breaks = breaks)
        }
        if (!is.na(xlim_lower) | !is.na(xlim_upper)) {
            figure <- figure + xlim(xlim_lower, xlim_upper)
        }
    }


    ##### Return #####
    model <- as.matrix(data.frame(q = q, xpos = "user",
                                  ci_type = ci, nboot = nboot))
    coefficients <- as.matrix(data.frame(beta = beta, se = se,
                                         lci = lci, uci = uci, pval = pval))
    rownames(coefficients) <- powers
    lace <- as.matrix(data.frame(beta = (frac_coef/xcoef),
                          se = (abs(frac_se/xcoef)),
                          lci = (frac_coef/xcoef - 1.96 * (abs(frac_se/xcoef))),
                          uci = (frac_coef/xcoef + 1.96 * (abs(frac_se/xcoef))),
                          pval = (2 * stats::rnorm(-abs(frac_coef / frac_se)))))
    rownames(lace) <- 1:nrow(lace)
    xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
    rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
    p_tests <- as.matrix(data.frame(fp_d1_d2 = p_d1_d2, fp = p_fp,
                                    quad = p_quadratic, Q = p_q))
    p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))

    if (fig == TRUE) {
        results <- list(n = NA, model = model, powers = powers,
                        coefficients = coefficients, lace = lace,
                        xcoef = xcoef_quant, p_tests = p_tests,
            p_heterogeneity = p_heterogeneity, figure = figure)
    }
    if (fig == FALSE) {
        results <- list(n = NA, model = model, powers = powers,
                        coefficients = coefficients, lace = lace,
                        xcoef = xcoef_quant, p_tests = p_tests,
            p_heterogeneity = p_heterogeneity)
    }
    class(results) <- "frac_poly_mr"
    return(results)
}
