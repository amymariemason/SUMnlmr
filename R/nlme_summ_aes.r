#Several of the functions here are taken from James Staley's nlmr package, to avoid
#need to install that package in order to create summarised data
# https://github.com/jrs95/nlmr for original functions

#' Hamardman product
#'
#' hamardman.prod computes the Hamardman product of a vector of regression
#' coefficients and a matrix of covariates.
#' @param coef vector of regression coefficients.
#' @param covar a matrix of covariates
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
hamardman.prod <- function(coef, covar){
  if(length(coef)!=ncol(covar)) stop("the number of coefficients is greater than
                                     the number of covariates")
  results <- reprow(coef,nrow(covar))*covar
  return(results)
}

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
#' @author Amy Mason <am2609@medschl.cam.ac.uk> based on similar function in nlmr
#' by James R Staley
iv_free <- function(y,
                    x,
                    g,
                    covar=NULL,
                    q=10,
                    family="gaussian",
                    controlsonly=T){
  if(family=="gaussian"){
    if(!is.null(covar)){
      model <- lm(x~g+covar)
    }else{
        model <- lm(x~g)
    }
    if(any(is.na(model$coef))) stop("there are missing regression coefficients
                                    in the regression of the exposure on the
                                    instrument and covariates")
    x0 <- resid(model)
  }
  if(family=="binomial"){
    if (controlsonly==T){
      if(!is.null(covar)){
        model <- lm(x[y == 0] ~ g[y == 0] + covar[y == 0, , drop = F])
        if(any(is.na(model$coef))) stop("there are missing regression coefficients
                                    in the regression of the exposure on the
                                    instrument and covariates in the controls")
        x0 <- x - (model$coef[1] + model$coef[2]*g +
                rowSums(hamardman.prod(model$coef[3:length(model$coef)],covar)))
      }else{
        model <- lm(x[y == 0] ~ g[y == 0])
        if(any(is.na(model$coef))) stop("there are missing regression
                                          coefficients in the regression of the
                                          exposure on the instrument and
                                        covariates in the controls")
      x0 <- x - (model$coef[1] + model$coef[2]*g)
      }
    }else{
      if(!is.null(covar)){
        model <- lm(x ~ g + covar)
        if(any(is.na(model$coef))) stop("there are missing regression coefficients
                                    in the regression of the exposure on the
                                    instrument and covariates in the controls")
        x0 <- x - (model$coef[1] + model$coef[2]*g +
                   rowSums(hamardman.prod(model$coef[3:length(model$coef)],covar)))
      }else{
        model <- lm(x ~ g)
        if(any(is.na(model$coef))) stop("there are missing regression
                                          coefficients in the regression of the
                                          exposure on the instrument and
                                        covariates in the controls")
        x0 <- x - (model$coef[1] + model$coef[2]*g)
      }
    }
  }
  xcoef <- model$coef[2]
  quantiles <- quantile(x0, probs=seq(0,1,1/q))
  x0q <- cut(x0, quantiles, include.lowest=T, labels=F)
  results <- list(xcoef=xcoef, x0=x0, x0q=x0q)
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
  p_heterogeneity <- as.data.frame(object$p_heterogeneity)
  if (is.null(object$figure)) {
    summ <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients, p_tests = p_tests,
      p_heterogeneity = p_heterogeneity
    )
  }
  if (!is.null(object$figure)) {
    summ <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients, p_tests = p_tests,
      p_heterogeneity = p_heterogeneity,
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
  cat("\n\nHeterogeneity tests")
  cat("\nCochran Q p-value:", signif(x$p_heterogeneity$Q, digits = 3))
  cat("\nTrend p-value:", signif(x$p_heterogeneity$trend, digits = 3))
  cat("\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}
