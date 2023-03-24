#' Instrumental variable analysis using piecewise linear method based on summary
#' data
#'
#' @description piecewise_summ_mr performs instumental variable analysis by
#' fitting piecewise linear functions to localised average causal effects
#'
#' @note The non-linearity tests uses 'method="DL"' to calculate the p-value for
#' the hetrogeneity trend. The fractional polynomial equivalent function allows
#' you to set the method, meaning you may get different results.
#' @note There is no option for covariates; they would need to be applied at an
#' earlier stage in the individual data, using the mr_summarise function.
#'
#' @param by vector of gene-outcome associations.
#' @param bx vector of gene-exposure associations.
#' @param byse vector of standard errors of gene-outcome associations.
#' @param bxse vector of standard errors of gene-exposure associations.
#' @param xmean average value of the original exposure in each iv-free strata
#' (or whatever summary of the exposure level in the stratum is desired).
#' @param xmin min value of the original exposure in each stratum (see note)
#' @param xmax max value of the original exposure in each stratum (see note)
#' @param xbreaks break points for the stratum x values (see note)
#' @param family a description of the error distribution and link function to be
#'  used in the model. This is a character string naming
#'  either the gaussian (i.e. "gaussian" for continuous data) or binomial
#'  (i.e. "binomial" for binary data) family function.
#' @param average.exposure.associations TRUE means that the bx estimates are averaged across strata, FALSE means that they are not. Default option is FALSE.
#' @param ci the type of 95\\% confidence interval. There are four options:
#' (i) using the model standard errors ('model_se'), (ii) using bootstrap
#' standard errors ('bootstrap_se'), (iii) using bootstrap percentile
#' confidence intervals ('bootstrap_per')
#' The default is the model standard errors.
#' @param nboot the number of bootstrap replications (if required). The default
#' is 1000 replications.
#' @param fig a logical statement as to whether the user wants the results
#' displayed in a figure. The default is false.
#' @param ref the reference point for the figure. This is the value of the
#' function that represents the expected difference in the outcome compared with
#'  this reference value when the exposure is set to different values. The
#'  default is the mean of x.
#' @param pref_x the prefix/label for the x-axis. The default is "x".
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis.
#'  The default is "x".
#' @param pref_y the prefix/label for the y-axis. The default is "y".
#' @param breaks breaks on the y-axis of the figure.
#' @param ci_fig setting confidence interval type. "point" places error bars
#' at the mean of each stratum; "line" draws upper and lower piecewise lines.
#' @param seed The random seed to use when generating the bootstrap samples (for reproducibility). If set to \code{NA}, the random seed will not be set.
#' @note The min and max of x stratum values are used to choose the appropiete
#' range for fitting of each causal estimate. In the code for summarising data,
#' this is set at the 10% point, and the 90% of each stratum; 20% and 80% in the
#'  external ends of the end strata. The first lower value and all upper value
#'  are used to set the break points for the estimates in the graph.
#'  Alternatively you can hardset this using xbreaks.
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
#' @author Amy Mason, leaning heavily on work by James Statley and Matt Arnold
#' @import ggplot2
#' @import stats
#' @importFrom matrixStats rowQuantiles
#' @importFrom metafor rma
#' @importFrom metafor rma.uni
#' @export
piecewise_summ_mr <- function(by,
                              bx,
                              byse,
                              bxse,
                              xmean,
                              xmin,
                              xmax,
                              xbreaks = NULL,
                              family = "gaussian",
                              average.exposure.associations = FALSE,
                              ci = "model_se",
                              nboot = 1000,
                              fig = T,
                              ref = mean(xmean),
                              pref_x = "x",
                              pref_x_ref = "x",
                              pref_y = "y",
                              breaks = NULL,
                              ci_fig = "point", seed=875) {

  if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
}
if (!is.na(seed)) { set.seed(seed) }


  ##### Error messages #####

  stopifnot(is.vector(by),
    is.vector(bx),
    is.vector(byse),
    is.vector(bxse),
    "by is not numeric" = (is.numeric(by) | is.integer(by)),
    "bx is not numeric" = (is.numeric(bx) | is.integer(bx)),
    "byse is not numeric" = (is.numeric(byse) | is.integer(byse)),
    "bxse is not numeric" = (is.numeric(bxse) | is.integer(bxse)),
    "xmean is not numeric" = (is.numeric(xmean) | is.integer(xmean)),
    "xmin is not numeric" = (is.numeric(xmin) | is.integer(xmin)),
    "xmax is not numeric" = (is.numeric(xmax) | is.integer(xmax)),
    "by is single value - cannot calculate multiple quantiles" =
      length(by) > 1,
    "difference number of observations for the outcome & exposure" =
      length(by) == length(bx),
    "missing by values" = !anyNA(by),
    "missing bx values" = !anyNA(bx),
    "missing byse values" = !anyNA(byse),
    "missing bxse values" = !anyNA(bxse),
    "missing xmean values" = !anyNA(xmean),
    "missing xmin values" = !anyNA(xmin),
    "missing xmax values" = !anyNA(xmax),
    "family has to be either gaussian or binomial" =
      family %in% c("gaussian", "binomial")
  )

  ##### define variables #######
  frac_coef <- by
  frac_se <- byse
  xcoef_sub <- bx
  xcoef_sub_se <- bxse
    if (average.exposure.associations == TRUE) { xcoef <- sum(bx * (bxse^(-2))) / sum(bxse^(-2)) }
   else { xcoef <- bx }
  q <- length(by)
  coef <- frac_coef / xcoef
  coef_se <- abs(frac_se / xcoef)

  ##### Non-linearity tests #####
  p_quadratic <- rma(coef ~ xmean, (coef_se)^2,
    method = "FE"
  )$pval[2]
  p_q <- 1 - pchisq(rma(coef, vi = (coef_se)^2)$QE, df = (q - 1))

  ##### Confidence Inteval ####
  if (ci == "bootstrap_per" | ci == "bootstrap_se") {
    boot_coef <- data.frame(matrix(ncol = q, nrow = nboot))
    for (i in 1:nboot) {
      # vary the value of by slightly
      boot_by <- by + rnorm(q, 0, byse)
      boot_bx <- bx + rnorm(q, 0, bxse)
      # recalculate the causal est based on this
      boot_xcoef <- sum(boot_bx * (bxse^(-2))) / sum(bxse^(-2))
      boot_coef[i, ] <- boot_by / boot_xcoef
    }
  }

  if (ci == "bootstrap_se") {
    cov <- var(boot_coef)
    se <- sqrt(diag(cov))
    lci <- coef - 1.96 * se
    uci <- coef + 1.96 * se
    pval <- 2 * pnorm(-abs(coef / se))
  }
  if (ci == "bootstrap_per") {
    se <- NA
    lci <- apply(boot_coef,
      MARGIN = 2,
      function(x) quantile(x, probs = 0.025)
    )
    uci <- apply(boot_coef,
      MARGIN = 2,
      function(x) quantile(x, probs = 0.975)
    )
    pval <- NA
  }
  if (ci == "model_se") {
    nboot <- NA
    se <- coef_se
    lci <- coef - 1.96 * coef_se
    uci <- coef + 1.96 * coef_se
    pval <- 2 * pnorm(-abs(coef / coef_se))
  }


  # estimate range creation

  if (is.null(xbreaks)) {
    quantiles<-NULL
    quantiles[1]<-xmin[1]
    quantiles[q+1]<-xmax[q]
    for (j in 2:q){
      quantiles[j]<-(xmin[j]+xmax[j-1])/2
    }
  } else {
    quantiles <- xbreaks
  }

  ##### Results #####
  lci <- as.numeric(lci)
  uci <- as.numeric(uci)


  ##### Figure #####
  if (fig == T) {
    figure <- piecewise_summ_figure(
      xcoef = xcoef,
      coef = coef,
      lci = lci,
      uci = uci,
      xmean = xmean,
      xbreaks = quantiles,
      ref = ref,
      pref_x = pref_x,
      family = family,
      pref_x_ref = pref_x_ref,
      pref_y = pref_y,
      breaks = breaks,
      ci_fig = ci_fig
    )
  }

  ##### Return #####
  model <- as.matrix(data.frame(q = q, nboot = nboot))
  lace <- as.matrix(data.frame(
    beta = coef,
    se = coef_se,
    lci = lci,
    uci = uci,
    pval = pval
  ))
  rownames(lace) <- seq_len(nrow(lace))
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- seq_len(nrow(xcoef_quant))
  p_tests <- as.matrix(data.frame(quad = p_quadratic, Q = p_q))

  if (fig == F) {
    results <- list(
      model = model,
      lace = lace,
      xcoef = xcoef_quant,
      p_tests = p_tests
    )
  } else {
    results <- list(
      model = model,
      lace = lace,
      xcoef = xcoef_quant,
      p_tests = p_tests,
      figure = figure
    )
  }
  class(results) <- "piecewise_summ_mr"
  return(results)
}


#' Piecewise linear figure
#'
#' piecewise_figure plots the piecewise linear function.
#' @import matrixStats
#' @import ggplot2
#' @param xcoef the association between the exposure and the instrument.
#' @param coef the coefficients of the localized causal effects.
#' @param xmean mean of each x stratum (for plotting point estimates)
#' @param lci upper confidence interval range for each causal effects
#' @param uci upper confidence interval range for each causal effects.
#' @param xbreaks breakpoints in x (for plotting line estimates)
#' @param ref the reference point for the figure. This is the value of the
#'  that represents the expected difference in the outcome compared with this
#'  reference value when the exposure is set to different values. The default
#'  is the mean of x.
#' @param pref_x the prefix/label for the x-axis. The default is "x".
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis.
#'  The default is "x".
#' @param pref_y the prefix/label for the y-axis. The default is "y".
#' @param breaks breaks on the y-axis of the figure.
#' @param family a description of the error distribution and link function to be
#'  used in the model. This is a character string naming
#'  either the gaussian (i.e. "gaussian" for continuous data) or binomial
#'  (i.e. "binomial" for binary data) family function.
#' @param ci_fig point confidence intervals, or as ribbon ("point" or "ribbon")
#' @return the plot of the piecewise linear function.
#' @author Amy Mason <am2609@medschl.cam.ac.uk>,
#' leaning on work by James Statley and Matt Arnold
#' @import ggplot2
#' @export
piecewise_summ_figure <- function(xcoef, coef,
                                  xmean, lci, uci,
                                  xbreaks,
                                  family = "gaussian",
                                  ref = mean(xmean),
                                  pref_x = "x",
                                  pref_x_ref = "x",
                                  pref_y = "y",
                                  breaks = NULL,
                                  ci_fig = "point") {

  # set what was variable number of quartiles in main function
  # forced by data here
  # q is number of quartiles for both lines and ci
  q <- length(coef)
  # these are the x breakpoints in the linear fit
  m <- xbreaks

  # counting variable
  l <- q + 1

  # local variables to stop R CMD Check complaining
  y_lci <- NULL
  y_uci <- NULL

  # find which section ref point is in
  for (i in 1:q) {
    if (m[i] <= ref & m[(i + 1)] >= ref) {
      ref_pos <- i + 1
    }
  }

  ### create the y coordinates for these quantiles breaks, based on the
  ### beta estimates in each quartile

  y_mm <- NULL
  uci_mm <- NULL
  lci_mm <- NULL

  y_mm[1] <- 0
  uci_mm[1] <- 0
  lci_mm[1] <- 0

  for (k in 2:l) {
    y_mm[k] <- (coef[k - 1] * m[k] - coef[k - 1] * m[k - 1]) + y_mm[k - 1]
    uci_mm[k] <- uci[k - 1] * m[k] - uci[k - 1] * m[k - 1] + uci_mm[k - 1]
    lci_mm[k] <- lci[k - 1] * m[k] - lci[k - 1] * m[k - 1] + lci_mm[k - 1]
  }
  # create reference points
  y_ref <- coef[ref_pos - 1] * ref -
    coef[ref_pos - 1] * m[ref_pos - 1] +
    y_mm[ref_pos - 1]

  uci_ref <- uci[ref_pos - 1] * ref -
    uci[ref_pos - 1] * m[ref_pos - 1] +
    uci_mm[ref_pos - 1]



  lci_ref <- lci[ref_pos - 1] * ref -
    lci[ref_pos - 1] * m[ref_pos - 1] +
    lci_mm[ref_pos - 1]

  # set ref points to zero
  y_mm_ref <- y_mm - y_ref

  uci_mm_ref <- pmax(uci_mm - uci_ref, lci_mm - lci_ref)
  lci_mm_ref <- pmin(uci_mm - uci_ref, lci_mm - lci_ref)

  # create y-cordinates for the mean of each segment
  y_mm_quant <- NULL
  lci_mm_quant <- NULL
  uci_mm_quant <- NULL

  for (j in 1:q) {
    x_ci <- xmean[j]
    # find segment containing the mean of the jth stratum
    for (i in 1:q) {
      if (m[i] <= x_ci & m[(i + 1)] >= x_ci) {
        ci_pos <- i + 1
      }
    }

    y_mm_quant[j] <- coef[ci_pos - 1] * x_ci -
      coef[ci_pos - 1] * m[ci_pos - 1] +
      y_mm[ci_pos - 1]
    lci_mm_quant[j] <- lci[ci_pos - 1] * x_ci -
      lci[ci_pos - 1] * m[ci_pos - 1] +
      lci_mm[ci_pos - 1]
    uci_mm_quant[j] <- uci[ci_pos - 1] * x_ci -
      uci[ci_pos - 1] * m[ci_pos - 1] +
      uci_mm[ci_pos - 1]
  }


  # rescale to ref point
  y_mm_quant_ref <- y_mm_quant - y_ref
  lci_mm_quant_ref <- pmin(lci_mm_quant - lci_ref,uci_mm_quant - uci_ref)
  uci_mm_quant_ref <- pmax(lci_mm_quant - lci_ref,uci_mm_quant - uci_ref)


  ##### Figure#####
  if (family != "binomial") {
    # collect data

    plot_data <- data.frame(
      x = m, y = y_mm_ref,
      y_lci = lci_mm_ref, y_uci = uci_mm_ref
    )
    plot_data1 <- data.frame(
      x = xmean, y = y_mm_quant_ref,
      y_lci = lci_mm_quant_ref,
      y_uci = uci_mm_quant_ref
    )
    plot_data2 <- data.frame(x = ref, y = 0)

    # set figure

    if (ci_fig == "point") {
      figure <- ggplot(plot_data, aes(x)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_line(aes(x = x, y = y), colour = "black") +
        geom_errorbar(aes(x = x, ymin = y_lci, ymax = y_uci),
          data = plot_data1,
          color = "grey", width = 0.025
        ) +
        geom_point(aes(x = x, y = y),
          data = plot_data1,
          colour = "black", size = 2
        ) +
        geom_point(aes(x = x, y = y),
          data = plot_data2,
          colour = "red", size = 2
        )
    } else {
      figure <- ggplot(plot_data, aes(x)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_ribbon(aes(ymin = y_lci, ymax = y_uci), alpha = 0.15) +
        geom_line(aes(x = x, y = y), colour = "black")
    }

    figure <- figure + theme_bw() + labs(
      x = pref_x,
      y = bquote(.(pref_y) ~ " [" ~ .(pref_x_ref)["ref"] ~ " = " ~ .(round(ref, 2)) ~ "]")
    ) +
      theme(
        axis.title.x = element_text(vjust = 0.5, size = 20),
        axis.title.y = element_text(vjust = 0.5, size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (!is.null(breaks)) {
      figure <- figure + scale_y_continuous(breaks = breaks)
    }
  } else {
    # collect data
    pref_y <- paste0("Odds ratio of ", pref_y)
    plot_data <- data.frame(
      x = m, y = exp(y_mm_ref), y_lci = exp(lci_mm_ref),
      y_uci = exp(uci_mm_ref)
    )
    plot_data1 <- data.frame(
      x = xmean, y = exp(y_mm_quant_ref),
      y_lci = exp(lci_mm_quant_ref),
      y_uci = exp(uci_mm_quant_ref)
    )
    plot_data2 <- data.frame(x = ref, y = 1)


    if (ci_fig == "point") {
      figure <- ggplot(plot_data, aes(x)) +
        geom_hline(aes(yintercept = 1), colour = "grey") +
        geom_line(aes(x = x, y = y), colour = "black") +
        geom_errorbar(aes(x = x, ymin = y_lci, ymax = y_uci),
          data = plot_data1,
          color = "grey", width = 0.025
        ) +
        geom_point(aes(x = x, y = y),
          data = plot_data1, colour = "black",
          size = 2
        ) +
        geom_point(aes(x = x, y = y),
          data = plot_data2, colour = "red",
          size = 2
        )
    } else {
      figure <- ggplot(plot_data, aes(x)) +
        geom_hline(aes(yintercept = 1), colour = "grey") +
        geom_ribbon(aes(ymin = y_lci, ymax = y_uci), alpha = 0.15) +
        geom_line(aes(x = x, y = y), colour = "black")
    }

    figure <- figure + theme_bw() + labs(
      x = pref_x,
      y = bquote(.(pref_y) ~ " [" ~ .(pref_x_ref)["ref"] ~ " =" ~ .(round(ref, 2)) ~ "]")
    ) +
      theme(
        axis.title.x = element_text(vjust = 0.5, size = 20),
        axis.title.y = element_text(vjust = 0.5, size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (!is.null(breaks)) {
      figure <- figure + scale_y_continuous(breaks = breaks)
    }
    ybreaks<-ggplot_build(figure)$layout$panel_params[[1]]$y$breaks
    figure<-figure + scale_y_continuous(trans="log", breaks=ybreaks)


  }
  return(figure)
}



#' Summary of piecewise linear fits
#'
#' summary method for class "piecewise_summ_mr".
#' @param object an object of class "piecewise_summ_mr".
#' @param ... Arguments to be passed to or from other methods,
#' @author Amy Mason <am2609@medschl.cam.ac.uk>
#' @export
summary.piecewise_summ_mr <- function(object, ...) {
  model <- as.data.frame(object$model)
  coefficients <- as.data.frame(object$lace)
  p_tests <- as.data.frame(object$p_tests)
  if (is.null(object$figure)) {
    summ <- list(
      model = model,
      coefficients = coefficients,
      p_tests = p_tests
    )
  }
  if (!is.null(object$figure)) {
    summ <- list(
      model = model,
      coefficients = coefficients,
      p_tests = p_tests,
      figure = object$figure
    )
  }
  class(summ) <- "summary.piecewise_summ_mr"
  return(summ)
}


#' Print summary of piecewise linear fits
#'
#' print summary method for class "piecewise_mr".
#' @param x an object of class "piecewise_mr".
#' @param ... Arguments to be passed to or from other methods,
#' @author Amy Mason <am2609@medschl.cam.ac.uk>
#' @export
print.summary.piecewise_summ_mr <- function(x, ...) {
  cat("Call: piecewise_summ_mr")
  cat("\n Quantiles: ", as.character(x$model$q), "; Number of bootstrap
      replications: ", as.character(x$model$nboot), sep = "")
  cat("\n\nLACE:\n")
  names(x$coefficients) <- c(
    "Estimate", "Std. Error", "95%CI Lower",
    "95%CI Upper", "p.value"
  )
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nNon-linearity tests")
  cat("\nQuadratic p-value:", signif(x$p_tests$quad, digits = 3))
  cat("\nCochran Q p-value:", signif(x$p_tests$Q, digits = 3))
  cat("\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}
