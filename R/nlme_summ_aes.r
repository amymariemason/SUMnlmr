#' Repeat rows
#'
#' This function creates a matrix of a repeated vector where each row is the same.
#' @param x vector to be repeated
#' @param n number of repeats
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}



#' Create plot of Fractional Polynomial Fit
#'
#' summary method for class "frac_poly_mr".
#' @param x an object of class "frac_poly_mr".
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
#' summary method for class "frac_poly_mr".
#' @param object an object of class "frac_poly_mr".
#' @param ... additional arguments affecting the summary produced.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.frac_poly_mr <- function(object, ...) {
  model <- as.data.frame(object$model)
  powers <- object$powers
  n <- object$n
  coefficients <- as.data.frame(object$coefficients)
  if (model$ci_type == "bootstrap_per") {
    coefficients <- coefficients[, c(1 ,3, 4)]
    }
  p_tests <- as.data.frame(object$p_tests)
  p_heterogeneity <- as.data.frame(object$p_heterogeneity)
  if (is.null(object$figure)) {
    summ <- list(model = model, powers = powers, n = n,
                 coefficients = coefficients, p_tests = p_tests,
                 p_heterogeneity = p_heterogeneity)
    }
  if (!is.null(object$figure)) {
    summ <- list(model = model, powers = powers, n = n,
                 coefficients = coefficients, p_tests = p_tests,
                 p_heterogeneity = p_heterogeneity, figure = object$figure)
    }
  class(summ) <- "summary.frac_poly_mr"
  return(summ)
}


#' Print Summary Fractional Polynomial Fits
#'
#' print.summary method for class "frac_poly_mr".
#' @param x an object of class "frac_poly_mr".
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
        as.character(x$model$q), "; 95%CI: ", ci_type, sep = "")
    }
  if(ci_type != "Model based SEs") {
    cat("\n\nNumber of individuals: ", x$n, "; Quantiles: ",
        as.character(x$model$q), "; 95%CI: ", ci_type,
        "; Number of bootstrap replications: ", as.character(x$model$nboot),
        sep = "")
    }
  cat("\n\nPowers:", x$powers)
  cat("\n\nCoefficients:\n")
  if (ci_type == "Percentile bootstrap") {
    names(x$coefficients) <- c("Estimate", "95%CI Lower", "95%CI Upper")
    stats::printCoefmat(x$coefficients)
    }
  if (ci_type != "Percentile bootstrap") {
    names(x$coefficients) <- c("Estimate", "Std. Error", "95%CI Lower",
                               "95%CI Upper", "p.value")
    stats::printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
    }
  cat("\nNon-linearity tests")
  cat("\nFractional polynomial degree p-value:", signif(x$p_tests$fp_d1_d2,
                                                        digits = 3))
  cat("\nFractional polynomial non-linearity p-value:", signif(x$p_tests$fp,
                                                               digits = 3))
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
