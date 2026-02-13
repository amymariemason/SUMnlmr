#Several of the functions here are taken from James Staley's nlmr package, to avoid
#need to install that package in order to create summarised data
# https://github.com/jrs95/nlmr for original functions

#' IV-free exposure
#'
#' iv_free computes the IV-free exposure.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param q the number of quantiles the exposure distribution is to be split
#' into. The default is deciles (i.e. 10
#' quantiles).
#' @param covar a matrix of covariates.
#' @param family a description of the error distribution and link function to be
#'  used in the model (either "gaussian" or "binomial" can be specified).
#' @param controlsonly whether to estimate the gx association in all people,
#' or in controls only. This is set to TRUE by default, but has no effect if
#' family is gaussian.
#' @return \item{xcoef}{the association between the exposure and the instrument}
#' @return \item{x0}{the IV-free exposure.}
#' @return \item{x0q}{the quantiles of x0.}
#' @export
#' @author Amy Mason <am2609@medschl.cam.ac.uk> based on similar function in nlmr
#' by James R Staley
iv_free <- function(y,
                    x,
                    g,
                    covar=NULL,
                    gxe_covar=NULL,
                    gxe_interaction=NULL,
                    q=10,
                    controlsonly=T){
  ##### build dataframe for model

  n <- length(x)
  if (length(g) != n || length(y) != n) stop("y, x, and g must have the same length.")

  df <- data.frame(y = y, x = x, g = g)

  # Add covariates (matrix/data.frame with columns)
  if (!is.null(covar)) {
    covar <- as.data.frame(covar)
    if (nrow(covar) != n) stop("covar must have the same number of rows as x.")
    # Avoid name clashes
    names(covar) <- make.names(names(covar), unique = TRUE)
    df <- cbind(df, covar)
  }

  # Add gxe_covariates (matrix/data.frame with columns)
  if (!is.null(gxe_covar)) {
    gxe_covar <- as.data.frame(gxe_covar)
    if (nrow(gxe_covar) != n) stop("covar must have the same number of rows as x.")
    # Avoid name clashes
    names(gxe_covar) <- make.names(names(gxe_covar), unique = TRUE)
    # prefix to avoid clashes with covars
    names(gxe_covar) <- paste0("cov_", names(gxe_covar))
    df <- cbind(df, gxe_covar)
  }

  # Add interaction variables (vector or matrix/data.frame)
  if (!is.null(gxe_interaction)) {
    gxe_interaction <- as.data.frame(gxe_interaction)
    if (nrow(gxe_interaction) != n) {
      stop("interaction must have the same number of rows as x.")
    }
    names(gxe_interaction) <- make.names(names(gxe_interaction), unique = TRUE)
    # prefix to avoid clashes with covars
    names(gxe_interaction) <- paste0("int_", names(gxe_interaction))
    df <- cbind(df, gxe_interaction)
  }

  ###### choose subset of data frame, if controls only selected

  fit_idx <- rep(TRUE, n)
  if (isTRUE(controlsonly)) {
    fit_idx <- (y == 0)
    if (!any(fit_idx)) stop("controlsonly=TRUE but there are no controls (y==0).")
  }

  #### build formula

  rhs_terms <- c("g")

  if (!is.null(covar) && !is.null(gxe_interaction)) {
    stop("cannot supply both covar and interaction matrices to iv_free")
  }

  if (!is.null(covar)) {
    rhs_terms <- c(rhs_terms, names(covar))
  }

  if (!is.null(gxe_covar)) {
    rhs_terms <- c(rhs_terms, names(gxe_covar))
  }

  if (!is.null(gxe_interaction)) {
    int_names <- names(gxe_interaction)
    # main effects + interactions with g
    rhs_terms <- c(rhs_terms, int_names, paste0("g:", int_names))
  }

  fml <- reformulate(rhs_terms, response = "x")

  # Note: how to pass error message if model fails to converge because same terms
  # passed in gxe_covar and gxe_interaction

  ### fit linear model

  model <- lm(fml, data = df, subset = fit_idx)

  if (anyNA(coef(model))) {
    stop("There are missing regression coefficients in the regression of the exposure on the instrument and covariates/interactions.")
  }

  ### predict residuals

  pred_all <- predict(model, newdata = df)
  x0 <- df$x - pred_all

  quantiles <- quantile(x0, probs=seq(0,1,1/q))
  # Ensure strictly increasing breaks in case many ties
  quantiles <- unique(quantiles)
  if (length(quantiles) < q+1) {
    stop("Residuals have too little variation to compute quantile bins.")
  }

  x0q <- cut(x0, quantiles, include.lowest=T, labels=F)
  results <- list(x0=x0, x0q=x0q)
  return(results)
}

#' Repeat rows
#'
#' This function creates a matrix of a repeated vector where each row is the
#' same.
#' @param x vector to be repeated
#' @param n number of repeats
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
reprow <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}



#' Create plot of Fractional Polynomial Fit
#'
#' summary method for class 'frac_poly_mr'.
#' @param x an object of class 'frac_poly_mr'.
#' @param ... additional arguments affecting the summary produced.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

print.frac_poly_mr <- function(x, ...) {
  cat("\nCall: \nfrac_poly_mr")
  cat("\n\nPowers:\n")
  cat(x$powers)
  cat("\n\nCoefficients:\n")
  cat(x$coefficients[, 1])
  cat("\n\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}



#' Summarizing Fractional Polynomial Fits
#'
#' summary method for class 'frac_poly_mr'.
#' @param object an object of class 'frac_poly_mr'.
#' @param ... additional arguments affecting the summary produced.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.frac_poly_mr <- function(object, ...) {
  model <- as.data.frame(object$model)
  powers <- object$powers
  n <- object$n
  coefficients <- as.data.frame(object$coefficients)
  if (model$ci_type == "bootstrap_per") {
    coefficients <- coefficients[, c(1, 3, 4)]
  }
  p_tests <- as.data.frame(object$p_tests)
  if (is.null(object$figure)) {
    summ <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients, p_tests = p_tests
    )
  }
  if (!is.null(object$figure)) {
    summ <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients, p_tests = p_tests,
      figure = object$figure
    )
  }
  class(summ) <- "summary.frac_poly_mr"
  return(summ)
}


#' Print Summary Fractional Polynomial Fits
#'
#' print.summary method for class 'frac_poly_mr'.
#' @param x an object of class 'frac_poly_mr'.
#' @param ... additional arguments affecting the summary produced.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.summary.frac_poly_mr <- function(x, ...) {
  cat("Call: frac_poly_mr")
  ci_type <- "Model based SEs"
  if (x$model$ci_type == "bootstrap_se") {
    ci_type <- "Bootstrap based SEs"
  }
  if (x$model$ci_type == "bootstrap_per") {
    ci_type <- "Percentile bootstrap"
  }
  if (ci_type == "Model based SEs") {
    cat("\n\nNumber of individuals: ", x$n, "; Quantiles: ",
      (x$model$q), "; 95%CI: ", ci_type,
      sep = ""
    )
  }
  if (ci_type != "Model based SEs") {
    cat("\n\nNumber of individuals: ", x$n, "; Quantiles: ",
      as.character(x$model$q), "; 95%CI: ", ci_type,
      "; Number of bootstrap replications: ",
      as.character(x$model$nboot),
      sep = ""
    )
  }
  cat("\n\nPowers:", x$powers)
  cat("\n\nCoefficients:\n")
  if (ci_type == "Percentile bootstrap") {
    names(x$coefficients) <- c("Estimate", "95%CI Lower", "95%CI Upper")
    stats::printCoefmat(x$coefficients)
  }
  if (ci_type != "Percentile bootstrap") {
    names(x$coefficients) <- c(
      "Estimate", "Std. Error",
      "95%CI Lower", "95%CI Upper", "p.value"
    )
    stats::printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
  }
  cat("\nNon-linearity tests")
  cat(
    "\nFractional polynomial degree p-value:",
    signif(x$p_tests$fp_d1_d2, digits = 3)
  )
  cat(
    "\nFractional polynomial non-linearity p-value:",
    signif(x$p_tests$fp, digits = 3)
  )
  cat("\nQuadratic p-value:", signif(x$p_tests$quad, digits = 3))
  cat("\nCochran Q p-value:", signif(x$p_tests$Q, digits = 3))
  cat("\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}
