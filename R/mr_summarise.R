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
#' This is a character string naming either the gaussian
#' (i.e. "gaussian" for continuous outcome data) or binomial (i.e. "binomial" for
#' binary outcome data) family function. "Coxph" can be used to fit survival data
#' - in this case y must be a Surv object.
#' @param controlsonly whether to estimate the gx association in all people,
#' or in controls only. This is set to FALSE as default.
#' It has no effect if family is set to "gaussian"
#' @param q the number of quantiles the exposure distribution is to be split
#' into. Within each quantile a causal effect will be fitted, known as a
#' localised average causal effect (LACE). The default is deciles (i.e. 10
#' quantiles).
#' @param prestrat the proportional size of pre-strata in the doubly-ranked
#' method. If prestrat = 1 (default), then pre-strata will contain
#' the number of individuals equal to the number of strata, and 1 individual
#' from each pre-stratum is selected into each stratum. If prestrat = 10, then
#' pre-strata contain 10 times the number of individuals as the number of
#' strata, and 10 individuals from each pre-stratum are selected into each
#' stratum. Larger pre-strata can improve the differentiation between
#' pre-strata, although if pre-strata are too large such that the instrument
#' values vary strongly within pre-strata, then the benefit of the doubly-ranked
#' method is lost.
#' @param strata_method what method to use for determining strata. By default
#' this is set to "ranked", using Haodong Tian's double ranked version to calculate
#' strata. The alternative is "residual" for determining the strata from the residual of the
#' exposure regressed on the instrument (As in Statley and Burgess paper). The
#' residual method relies on a constant relationship between the instrument and the
#' exposure across the range of the exposure.
#' @param strata_bound controls what range to use for the LACE estimates in graphs display.
#' By default this is set to restricted, taking the 10th and 90th percentile of
#' internal strata and the 20th and 80th for the bottom of the lowest strata and
#' top of the highest strata. It is a vector taking the percentiles for the
#' lowers bounds of the bottom and then other strata and then upper bounds of
#' top and other strata.
#' This only impacts the "max" and "min" values for the summary table
#' This can be overridden in piecewise_summ_mr by using the xbreaks argument to
#' hardset different breakpoints or replacing default with c(0,0,1,1) to return
#' to true max and minimum
#' @param extra_statistics This will add a second output reporting extra
#' statistics for each strata. These include the true max and min of each
#' strata (regardless of strata_bound setting) and the f statistic and p-value
#' for the regressions
#' @param report_GR This will add the Gelman-Rubin statistics for each strata
#' to the output. Note this only works if strata_method="ranked".
#' @param report_het This will add p-values for assessing the heterogeneity of
#' the instrument - exposure relationship.
#' The first column is the p-value of the Cochran Q heterogeneity test (Q);
#' the second column is the p-value from the trend test (trend).
#' @param seed The random seed to use when generating the quantiles (for reproducibility). If set to \code{NA}, the random seed will not be set.
#' @return model the model specifications. The first column is the number of
#' quantiles (q); the second column is the position used to relate x to the LACE
#'  in each quantiles (xpos); the third column is the type of confidence
#'  interval constructed (ci); the fourth column is the number of bootstrap
#'  replications performed (nboot).

#' @author Amy Mason, leaning heavily on work by James Statley and Matt Arnold
#' @import ggplot2
#' @import matrixStats
#' @importFrom dplyr mutate group_by row_number arrange
#' @importFrom stats quantile
#' @importFrom metafor rma
#' @importFrom metafor rma.uni
#' @importFrom survival coxph
#' @importFrom survival is.Surv
#' @export
create_nlmr_summary <- function(y,
                                x,
                                g,
                                covar = NULL,
                                family = "gaussian",
                                controlsonly=FALSE,
                                q,
                                prestrat=1,
                                strata_method="ranked",
                                strata_bound=c(0.2,0.1,0.8,0.9),
                                extra_statistics =FALSE,
                                report_GR=FALSE,
                                report_het=FALSE,
                                seed=1234) {

  if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
}
if (!is.na(seed)) { set.seed(seed) }


  # checks
  stopifnot(
    "report_GR only works with strata_method ranked" = !(report_GR==TRUE &
                                                      strata_method!="ranked")
  )
  # coxph checks
  stopifnot(
    "y must be a Surv object with family coxph" = !(family=="coxph" &
                                                           !is.Surv(y))
  )
  stopifnot(
    "cannot use controlsonly option with family coxph" = !(family=="coxph" &
                                                      controlsonly==TRUE)
  )

  # covar issue
  if (!is.null(covar) & !(is.matrix(covar) & is.numeric(covar))) {
    warning("covariates should be entered as numeric matrix:
            attempting to covert", immediate.=TRUE)
    covar2<-model.matrix(~.,data=as.data.frame(covar), na.action=na.pass)[,-1]
    print(head(covar2))
    user_input <- readline("Do you want to run using this matrix for covariates (y/n) ")
    if(user_input != 'y') stop('Exiting since you did not press y')
    covar<-as.matrix(covar2)

  }


  # calculate the iv-free association
  if (family=="binomial" |family=="gaussian"| family=="coxph") {
  if (strata_method=="residual"){
    family2= ifelse(family=="binomial", "binomial", "gaussian")
    ivf <- iv_free(
    y = y, x = x, g = g,
    covar = covar, q = q, family = family2, controlsonly=controlsonly
  )
    x0q <- ivf$x0q
  } else if(strata_method=="ranked") {
    # haodong ranked strata method
    z = rank(g, ties.method = "random")
    strata1 = floor((z-1)/q/prestrat)+1
    # check GR statistic
    GR_stats<-getGRvalues(X=x, Zstratum=strata1)

    id= seq(x)
    temp<- data.frame(x=x,strata1=strata1,id=id, g=g)
    temp<- arrange(.data=temp, x)
    temp<-group_by(.data=temp, strata1)
    temp<-mutate(.data=temp, x0q= ceiling(rank(x, ties.method = "random")/prestrat))
    temp<-arrange(.data=temp, id)

    x0q <- temp$x0q
  } else {
    stop("strata ordering must be ranked or residual")
  }
  } else {
    stop("family must be gaussian or binomial or coxph")
  }


  quant <- q

  # this calculates the association for each quanta
  by <- rep(NA, quant)
  byse <- rep(NA, quant)
  bx <- rep(NA, quant)
  bxse <- rep(NA, quant)
  xmean <- rep(NA, quant)
  xmax <- rep(NA, quant)
  xmin <- rep(NA, quant)
  true_xmax <- rep(NA, quant)
  true_xmin <- rep(NA, quant)
  strata_stats<-vector("list", quant)

#use the ivfree quantiles

  for (j in 1:quant) {
    # describe the quantiles of original data
    if (j==1){
      xmin[j] <- quantile(x[x0q == j], strata_bound[1])
    }else{
      xmin[j] <- quantile(x[x0q == j], strata_bound[2])
    }
    if (j==quant){
      xmax[j] <- quantile(x[x0q == j], strata_bound[3])
    }else{
      xmax[j] <- quantile(x[x0q == j], strata_bound[4])
    }
    xmean[j] <- mean(x[x0q == j])
    # model the beta coefficients
    if (family == "binomial") {
      if (is.null(covar)) {
        model <- glm(y[x0q == j] ~ g[x0q == j], family = "binomial")
        if (controlsonly==T){
          model2 <- lm(x[x0q == j& y==0] ~ g[x0q == j & y==0])
        } else {
          model2 <- lm(x[x0q == j] ~ g[x0q == j])
        }
      }else{
        model <- glm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F],
                     family = "binomial")
        if (controlsonly==T){
          model2 <- lm(x[x0q == j& y==0] ~ g[x0q == j & y==0]+ covar[x0q == j & y==0, , drop = F])
        } else {
          model2 <- lm(x[x0q == j] ~ g[x0q == j, drop = F] + covar[x0q == j, , drop = F])
          }
      }
      if (is.na(model$coef[2])) {
        stop("the regression coefficient of the outcome on the instrument
           in one of the quantiles is missing")
      }
      by[j] <- model$coef[2]
      byse[j] <- summary(model)$coef[2, 2]
    }else if(family == "coxph"){
      if (is.null(covar)) {
        model <- coxph(y[x0q == j] ~ g[x0q == j])
        model2 <- lm(x[x0q == j] ~ g[x0q == j])
      }else{
        model <- coxph(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F])
        model2 <- lm(x[x0q == j] ~ g[x0q == j, drop = F] + covar[x0q == j, , drop = F])
      }
      if (is.na(model$coef[1])) {
        stop("the regression coefficient of the outcome on the instrument
           in one of the quantiles is missing")
      }
      by[j] <- model$coef[1]
      byse[j] <- summary(model)$coef[1, 3]

    }else {
      if (is.null(covar)) {
        model <- lm(y[x0q == j] ~ g[x0q == j])
        model2 <- lm(x[x0q == j] ~ g[x0q == j])
      }else{
        model <- lm(y[x0q == j] ~ g[x0q == j]+   covar[x0q == j, , drop = F])
        model2 <- lm(x[x0q == j] ~ g[x0q == j]+   covar[x0q == j, , drop = F])
      }
      if (is.na(model$coef[2])) {
        stop("the regression coefficient of the outcome on the instrument
           in one of the quantiles is missing")
      }
      by[j] <- model$coef[2]
      byse[j] <- summary(model)$coef[2, 2]
    }

if (is.na(model2$coef[2])) {
      stop("the regression coefficient of the exposure on the instrument
           in one of the quantiles is missing")
    }

       bx[j] <- model2$coef[2]
    bxse[j] <- summary(model2)$coef[2, 2]


    if (extra_statistics) {
      stats<- list( strata=j,
        xmin = quantile(x[x0q == j], 0),
                         xmax = quantile(x[x0q == j], 1),
                   xmean = mean(x[x0q == j]),
                         ymin = quantile(y[x0q == j], 0),
                         ymax = quantile(y[x0q == j], 1),
                         x_fstat = (bx[j]/bxse[j])^2

                   )
      strata_stats[[j]]<-append(strata_stats[[j]], stats)

    }


    model <- NULL
    model2 <- NULL
  }
  # output data

  output <- data.frame(bx, by, bxse, byse, xmean, xmin, xmax)
  names(output) <- c("bx", "by", "bxse", "byse", "xmean", "xmin", "xmax")

  final_output_list= list(summary=output)
  if (extra_statistics) {
    stats<- as.data.frame(do.call(rbind, strata_stats))
    final_output_list[["strata_statistics"]]<- stats
  }
  if (strata_method=="ranked"){
    final_output_list[["GR_max"]]<- GR_stats[1]

    ##### Test of IV-exposure assumption #####
    xcoef_sub <- bx
    xcoef_sub_se <- bxse
    p_het <- 1 - pchisq(rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE,
                        df = (q - 1)
    )
    p_het_trend <- rma.uni(xcoef_sub ~ xmean,
                           vi = xcoef_sub_se^2,
                           method = "DL"
    )$pval[2]

  if (report_GR==TRUE){
    final_output_list[["GR_results"]]<-GR_stats
  }
  if (report_het==TRUE){
    p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
    final_output_list[["Heterogeneity_results"]]<- p_heterogeneity
  }
    }


  # print(list(summary = head(output)))
    invisible(final_output_list)

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
#' @param family a description of the error distribution and link function to be used in the model.
#' For piecewise_mr this can be a character string naming either the gaussian
#' (i.e. "gaussian" for continuous data) or binomial (i.e. "binomial" for
#' binary data) family function.
#' @param controlsonly whether to estimate the gx association in all people,
#' or in controls only. This is set to TRUE as default.
#' It has no effect if family is set to "gaussian"
#' @param q The number of quantiles.
#' @param keep Whether to retain the individual level data as well as the summary data
#' @return summary A data-frame containing the semi-summarised beta_X and
#' beta_Y values, mean of X_0, mean of X, max of X and min of X for each quantile.
#' @author Amy Mason
#' @import stats
#' @export

create_summary_data <- function(Ytype, q = 10, keep = FALSE,
                                family = "gaussian",
                                controlsonly="TRUE", ...) {
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
    family = family,
    controlsonly=controlsonly
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

