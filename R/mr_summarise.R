#' Creation of summarised mendelian randomisation local estimates
#'
#' @description create_nlmr_summary takes individual level data and creates summerised
#' dataset, ready to save and share for summarised nlmr
#'
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param covar a matrix of covariates.
#' @param family a description of the error distribution and link function to be used in the model.
#' For piecewise_mr this can be a character string naming either the gaussian
#' (i.e. "gaussian" for continuous data) or binomial (i.e. "binomial" for
#' binary data) family function.
#' @param q the number of quantiles the exposure distribution is to be split
#' into. Within each quantile a causal effect will be fitted, known as a
#' localised average causal effect (LACE). The default is deciles (i.e. 10
#' quantiles).
#' @return model the model specifications. The first column is the number of
#' quantiles (q); the second column is the position used to relate x to the LACE
#'  in each quantiles (xpos); the third column is the type of confidence
#'  interval constructed (ci); the fourth column is the number of bootstrap
#'  replications performed (nboot).

#' @author Amy Mason, leaning heavily on work by James Statley and Matt Arnold
#' @import ggplot2
#' @import matrixStats
#' @importFrom nlmr iv_free
#' @importFrom stats quantile
#' @export
create_nlmr_summary <- function(y,
                                x,
                                g,
                                covar = NULL,
                                family = "gaussian",
                                q) {

  # calculate the iv-free association
  ivf <- iv_free(
    y = y, x = x, g = g,
    covar = covar, q = q, family = family
  )

  x0q <- ivf$x0q
  quant <- q

  # this calculates the association for each quanta
  by <- rep(NA, quant)
  byse <- rep(NA, quant)
  bx <- rep(NA, quant)
  bxse <- rep(NA, quant)
  x0mean <- rep(NA, quant)
  xmean <- rep(NA, quant)
  xmax <- rep(NA, quant)
  xmin <- rep(NA, quant)

  # find upper and lower cutoffs for each quartile of x NOT IVFREE
  quantiles_x <- quantile(x, probs = seq(0, 1, 1 / quant))
  x_quantiles <- cut(x, quantiles_x, include.lowest = T)
  x_quantiles <- as.numeric(x_quantiles)

  for (j in 1:quant) {
    # describe the quantiles of original data
    xmin[j] <- quantile(x[x_quantiles == j], 0.000001)
    xmax[j] <- quantile(x[x_quantiles == j], 0.999999)
    xmean[j] <- mean(x[x_quantiles == j])
    # model the y coefficient
    if (family == "gaussian") {
      if (is.null(covar)) {
        model <- lm(y[x0q == j] ~ g[x0q == j])
      } else {
        model <- lm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F])
      }
    }
    if (family == "binomial") {
      if (is.null(covar)) {
        model <- glm(y[x0q == j] ~ g[x0q == j], family = "binomial")
      } else {
        model <- glm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F],
          family = "binomial"
        )
      }
    }
    if (is.na(model$coef[2])) {
      stop("the regression coefficient of the outcome on the instrument
           in one of the quantiles is missing")
    }
    by[j] <- model$coef[2]
    byse[j] <- summary(model)$coef[2, 2]
    model <- NULL
    # model the x coefficient
    if (is.null(covar)) {
      model2 <- lm(x[x0q == j] ~ g[x0q == j])
    } else {
      model2 <- lm(x[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F])
    }
    bx[j] <- model2$coef[2]
    bxse[j] <- summary(model2)$coef[2, 2]
    model2 <- NULL
    x0mean[j] <- mean(x[x0q == j])
  }
  # output data
  output <- data.frame(bx, by, bxse, byse, x0mean, xmean, xmin, xmax)
  names(output) <- c("bx", "by", "bxse", "byse", "x0mean", "xmean", "xmin", "xmax")
  # print(list(summary = head(output)))
  invisible(list(summary = output))
}

#' Generation of individual level data
#' @description generates individual level data with a single genetic variant
#'
#' @param N number of individuals to create
#' @param gpar genetic parameter; used to create g: single genetic snp, from a
#' binomial distribution with n=2 and p = gpar.
#' @param par1 power parameter for fractional poly generation. See details
#' @param par2 power parameter for fractional poly generation. . See details.
#' @param beta0 covariate parameter. See details
#' @param beta1 covariate parameter. See details
#' @param beta2 covariate parameter. See details
#' @param confound confounding parameter,c. See details.
#' @note This function generates a database with genetic relationships suitable
#' for evaluating non-linear MR relationships.
#' A unknown covariate,u,  is generated as a N(0,1) variable.
#' Error terms are generated: Ex ~exp(1) and for Ey ~ N(0,1)
#' \eqn{X= 2+ 0.25*g +u + E_x}
#' Outcomes are as follows
#' \itemize{
#' \item Linear: \eqn{Y=b_0+ b_1 X + cU +E_y}
#' \item Quadratic \eqn{Y = b_0 + b_1 X + b_2 X^2 + cU +E_y}
#' \item Squareroot \eqn{Y = b_0 + b_1 \sqrt{X} + cU +E_y}
#' \item Log \eqn{Y = b_0 + b_1 \log(X) + cU +E_y}
#' \item Threshold \eqn{Y = b_0+ b_1 X + cU +E_y} if\eqn{X>b_2} and
#' \eqn{Y = b_0 + cU +E_y} otherwise
#' \item fracpoly \eqn{Y = b_0 + b_1 X^{p_1} + b_2 X^{p_2} + cU + E_y}
#'  with the usual adaptions for p=0 or p_1=p_2
#'  }
#' @return data A data-frame containing the values of g, the genetic variate;
#' X, the exposure; and a variety of Y, the outcome values.
#' All outcomes are continuous not binary.
#' @author Amy Mason
#' @import stats
#' @export
create_ind_data <- function(N, gpar = 0.3, par1 = 1, par2 = 0,
                            beta0 = 0, beta1 = 3, beta2 = 7, confound = 0.8) {
  # generate G
  data <- as.data.frame(rbinom(N, 2, gpar))
  names(data) <- c("g")

  # generate Unknown confound
  data$u <- runif(N, 0, 1)

  # generate error terms
  data$errorX <- rexp(N, 1)
  data$errorY <- rnorm(N, 0, 1)

  # build X
  data$X <- 2 + 0.25 * data$g + data$u + data$errorX

  # generate various Y with different exposure-outcome results
  data$linear.Y <- beta0 + beta1 * data$X + confound * data$u + data$errorY
  data$quadratic.Y <- beta0 + beta2 * (data$X)^2 + beta1 * data$X + confound * data$u +
    data$errorY
  data$sqrt.Y <- beta0 + beta1 * sqrt(data$X) + confound * data$u + data$errorY
  data$log.Y <- beta0 + beta1 * log(data$X) + confound * data$u + data$errorY
  data$threshold.Y <- ifelse(data$X > beta2, beta0 + beta1 * data$X, beta0) +
    confound * data$u + data$errorY
  if (par1 == par2) {
    if (par1 == 0) {
      data$fracpoly.Y <- beta0 + beta1 * log(data$X) + beta2 * log(data$X) * log(data$X) +
        confound * data$u + data$errorY
    } else {
      data$fracpoly.Y <- beta0 + beta1 * data$X^par1 + beta2 * log(data$X) * data$X^par1 +
        confound * data$u + data$errorY
    }
  } else {
    if (par1 == 0) {
      data$fracpoly.Y <- beta0 + beta1 * log(data$X) + beta2 * data$X^par2 +
        confound * data$u + data$errorY
    } else if (par2 == 0) {
      data$fracpoly.Y <- beta0 + beta1 * data$X^par1 + beta2 * log(data$X) +
        confound * data$u + data$errorY
    } else {
      data$fracpoly.Y <- beta0 + beta1 * data$X^par1 + beta2 * data$X^par2 +
        confound * data$u + data$errorY
    }
  }
  return(data)
}


#' generation of summary level data
#' @description create_summary_data generates semi-summarized level data
#' @param Ytype The relationship between X and Y; can be "linear", "quad", "sqrt", "log", "threshold" or "fracpoly"
#' @param ... parameters passed to mr_create_data for control of X-Y relationship
#' @param q The number of quantiles.
#' @param keep Whether to retain the individual level data as well as the summary data
#' @return summary A data-frame containing the semi-summarised beta_X and
#' beta_Y values, mean of X_0, mean of X, max of X and min of X for each quantile.
#' @author Amy Mason
#' @import stats
#' @export

create_summary_data <- function(Ytype, q = 10, keep = FALSE, ...) {
  # create Ytype_name to generate the appropriate function type
  if (Ytype == "linear") {
    Ytype_name <- "linear.Y"
  } else if (Ytype == "quad") {
    Ytype_name <- "quadratic.Y"
  } else if (Ytype == "sqrt") {
    Ytype_name <- "sqrt.Y"
  } else if (Ytype == "log") {
    Ytype_name <- "log.Y"
  } else if (Ytype == "threshold") {
    Ytype_name <- "threshold.Y"
  } else if (Ytype == "fracpoly") {
    Ytype_name <- "fracpoly.Y"
  } else {
    stop("model type not supported")
  }

  # create the data

  data <- create_ind_data(...)
  y <- data[, Ytype_name]
  g <- data$g
  x <- data$X
  summ <- create_nlmr_summary(
    y = y,
    x = x,
    g = g,
    q = q,
    family = "Gaussian",
    ...
  )
  summ_data <- summ$summary

  # keep entire set if keep variable set to TRUE
  if (keep == TRUE) {
    invisible(list(
      summary = summ_data,
      alldata = data
    ))
  } else {
    invisible(list(summary = summ_data))
  }
}
