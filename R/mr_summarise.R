#' Creation of summarised mendelian randomisation local estimates
#'
#' @description create_nlmr_summary takes individual level data and creates summerised
#' dataset, ready to save and share for summarised nlmr
#'
#' @details
#' ## Stratification models
#'
#' Given vectors of an exposure `x`, an outcome `y` and a instrument `g`, as well as a
#' matrices of covariates `E`. This package implements
#' several possible models for stratifying over the exposure value.
#'
#' The simplest of these the residual method, as described in Statley and Burgess,
#' 2017 <doi:https://doi.org/10.1002/gepi.22041>.  where we regress `x` on the genetic
#' instrument and the covariate matrix `E`
#' \deqn{x = \beta_0 + \beta_g g + \beta_E E}
#' where `\beta_E` is a suitable vector of coefficients, and calculate the residual value
#' \deqn{x_{resid} = x- (\beta_0 + \beta_g G + \beta_E E)}
#' Once the strata are formed, the genetic associations with `x` and `y` are calculated,
#' including the covariates matrix E e.g.
#' \deqn{x = \beta_1 + \beta_x g + \beta_{E1} E}
#' \deqn{y = \beta_2 + \beta_y g + \beta_{E2} E}
#' with \eqn{\beta_x} and \eqn{\beta_y} returned for each strata.
#' This method performs poorly unless the genetic effects on the exposure are constant.
#' This assumption is often implausible, and thus this approach is not recommended.
#' See Burgess, 2023 <doi:https://doi.org/10.1159/000531659>.
#'
#' The second option is "ranked", using Haodong Tian's double ranked version to calculate
#' strata. See Haodong et al, 2022 <doi: https://doi.org/10.1101/2022.06.28.497930> for more details on this calculation.
#' However, there are concerns about this method's application, particularly in UK Biobank,
#' with thanks to Hamiliton et al, 2023 <doi: https://doi.org/10.1007/s10654-024-01113-9> for their
#' examples of this. In particular, this method performs poorly when there is a
#' GxE interaction, as explained in Zhao et al, 2026 <doi: https://doi.org/10.64898/2026.01.22.26344640>
#'
#'The final method "interaction" is an extension of the ranked method developed
#'specifically to mitigate GxE-induced bias in that method, as detailed in
#'Zhou et al, 2026 <doi: https://doi.org/10.64898/2026.01.22.2634464>.
#'In this approach, an interaction model is first fitted using two
#'additional matrices, passed using the `gxe_covar` (F, covariates) and
#'`gxe_interaction` (H, effect modifiers).
#' \deqn{x = \beta_0 + \beta_g g + \beta_{F} F + \beta_{H} H + \beta_{g \times H} (g \times H)}
#' Then an interaction corrected value of X, \eqn{X - \beta_(g \times H)(g \times H)}
#' is used to form strata using the ranked method.
#' Within the strata, the associations with exposure and outcome are then
#' calculated using the usual `covar` covariants matrix E e.g.
#' \deqn{x = \beta_1 + \beta_x g + \beta_{E1} E}
#' \deqn{y = \beta_2 + \beta_y g + \beta_{E2} E}
#' with \eqn{\beta_x} and \eqn{\beta_y} returned for each strata.
#'
#' The two matrices supplied in `gxe_covar` and `gxe_interaction` are used to allow the
#' greatest flexibility in choice of models for correcting potential GxE interactions.
#'
#'
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param covar an optional matrix of covariates used to derive the stratification
#' and the genetic associations. If `interaction` is also provided, then these
#' are only used in the genetic association calculation.
#' @param gxe_covar an optional matrix of covariates used to derive the
#' stratification only. See details.
#' @param gxe_interaction an optional matrix of covariates used to derive the
#' stratification only. See details.
#' @param strata_method what method to use for determining strata. There are
#' three options "residual", "ranked", "interaction".
#' See details for a longer explanation of these methods. By default
#' this is set to "ranked".
#' @param x_corrected an optional numeric vector used if method is "ranked". This will
#' override the model and use these values in the ranked method stratification
#' calculation. This can also be used as a way of creating strata with your own
#' interaction model - supply the interaction corrected X values here.
#' @param x_strata an optional numeric vector used to replace the stratification
#' calculation. Only use this if you have precalculated the strata and just want
#' the genetic associations within those strata.
#' @param q the number of quantiles the exposure distribution is to be split
#' into. Within each quantile a causal effect will be fitted, known as a
#' localised average causal effect (LACE). The default is deciles (i.e. 10
#' quantiles).
#' @param family a description of the error distribution and link function to be
#' used in the model. It must be one of "gaussian", "binomial" or "coxph".
#' Gaussian should be used for continuous outcome data and will fit linear
#' regression models. Binomial should be used to fit logistic regression models to
#' binary outcome data. Coxph should be used to fit will fit binary outcome
#' This is a character string naming either the "gaussian"
#' (i.e. "gaussian" for continuous outcome data) or binomial (i.e. "binomial" for
#' binary outcome data) family function. "Coxph" can be used to fit survival data
#' - in this case y must be a Surv object.
#' @param controlsonly Only applied if family is "binomial" or "coxph".
#' If true, the genetic association with x is only calculated in the controls only.
#' @param prestrat Only applied if method is "ranked".
#' The proportional size of pre-strata in the doubly-ranked method.
#' If prestrat = 1 (default), then pre-strata will contain
#' the number of individuals equal to the number of strata, and 1 individual
#' from each pre-stratum is selected into each stratum. If prestrat = 10, then
#' pre-strata contain 10 times the number of individuals as the number of
#' strata, and 10 individuals from each pre-stratum are selected into each
#' stratum. Larger pre-strata can improve the differentiation between
#' pre-strata, although if pre-strata are too large such that the instrument
#' values vary strongly within pre-strata, then the benefit of the doubly-ranked
#' method is lost.
#' @param strata_bound This controls what range to use for the LACE estimates in
#' graphs display. By default, this is taken conservatively with the 10th and
#' 90th percentile of internal strata and the 20th and 80th for the bottom of the
#' lowest strata and top of the highest strata. It is supplied as a vector of percentiles
#' (lower bottom strata, lower other, higher top strata, higher other).
#' This only impacts the "max" and "min" values for the summary table.
#' This can be overridden in `piecewise_summ_mr()` by using the `xbreaks` argument to
#' hardset different breakpoints or replacing default with `c(0,0,1,1)` to return
#' to true max and minimums.
#' @param extra_statistics This will add a second output reporting extra
#' statistics for each strata. These include the true max and min of each
#' strata (regardless of strata_bound setting) and the f statistic and p-value
#' for the regressions.
#' @param report_GR This will add the Gelman-Rubin statistics for each strata
#' to the output. Note this only works if the strata method is "ranked".
#' @param report_het This will add p-values for assessing the heterogeneity of
#' the instrument - exposure relationship.
#' The first column is the p-value of the Cochran Q heterogeneity test (Q);
#' the second column is the p-value from the trend test (trend).
#' @param return_assignment reports the strata assignment of each participant,
#' to assist in examining stratum distributions.
#' @param seed The random seed to use when generating the quantiles
#' (for reproducibility). If set to \code{NA}, the random seed will not be set.
#' @return A list containing summary data and optional diagnostics.
#' \itemize{
#' \item \code{summary}: data frame with one row per stratum containing
#' \code{bx}, \code{by}, \code{bxse}, \code{byse}, \code{xmean}, \code{xmin}, and
#' \code{xmax}.
#' \item \code{Dropped}: always returned; integer giving the number of
#' individuals removed due to missing values in any of \code{y}, \code{x},
#' \code{g}, \code{covar}, \code{gxe_covar}, \code{gxe_interaction},
#' \code{x_corrected}, or \code{x_strata}. A message is also printed when
#' this value is greater than zero.
#' \item \code{strata_statistics}: returned when
#' \code{extra_statistics = TRUE}; per-stratum min/max of \code{x} and \code{y},
#' mean \code{x}, and exposure-model F statistic.
#' \item \code{GR_max}: maximum Gelman--Rubin statistic for ranked strata.
#' \item \code{GR_results}: full Gelman--Rubin output when
#' \code{report_GR = TRUE} and \code{strata_method = "ranked" or "interaction"}.
#' \item \code{Heterogeneity_results}: heterogeneity/trend test p-values when
#' \code{report_het = TRUE} and \code{strata_method = "ranked" or "interaction"}.
#' \item \code{Strata_assignment}: strata assignment for each participant when
#' return_assignment = TRUE
#' }
#' @author Amy Mason
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
                                gxe_covar=NULL,
                                gxe_interaction=NULL,
                                strata_method="ranked",
                                x_corrected =NULL,
                                x_strata=NULL,
                                family = "gaussian",
                                controlsonly=FALSE,
                                q=10,
                                prestrat=1,
                                strata_bound=c(0.2,0.1,0.8,0.9),
                                extra_statistics =FALSE,
                                report_GR=FALSE,
                                report_het=FALSE,
                                return_assignment=FALSE,
                                seed=1234) {

  if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
}
if (!is.na(seed)) { set.seed(seed) }

##################### in put checks
  # checks on method choice
  stopifnot("strata method must be one of residual, ranked or interaction" =
              strata_method%in% c("residual", "ranked", "interaction"))

  stopifnot(
    "report_GR does not work with strata method residual" =
      !(report_GR==TRUE & strata_method=="residual")
  )
  stopifnot("family must be one of gaussian, binomial or coxph" = family%in% c("gaussian", "binomial", "coxph"))

    # covar entry issue
  if (!is.null(covar) & !(is.matrix(covar) & is.numeric(covar))) {
    warning("covariates should be entered as numeric matrix:
            attempting to covert", immediate.=TRUE)
    covar2<-model.matrix(~.,data=as.data.frame(covar), na.action=na.pass)[,-1]
    print(head(covar2))
    user_input <- readline("Do you want to run using this matrix for covariates (y/n) ")
    if(user_input != 'y') stop('Exiting since you did not press y')
    covar<-as.matrix(covar2)

  }
  # interaction entry issue
  if (!is.null(gxe_interaction) & !(is.matrix(gxe_interaction) & is.numeric(gxe_interaction))) {
    warning("covariates should be entered as numeric matrix:
            attempting to covert", immediate.=TRUE)
    interaction2<-model.matrix(~.,data=as.data.frame(gxe_interaction), na.action=na.pass)[,-1]
    print(head(interaction2))
    user_input <- readline("Do you want to run using this matrix for interaction terms (y/n) ")
    if(user_input != 'y') stop('Exiting since you did not press y')
    gxe_interaction<-as.matrix(interaction2)

  }

  # second covar entry issue
  if (!is.null(gxe_covar) & !(is.matrix(gxe_covar) & is.numeric(gxe_covar))) {
    warning("covariates should be entered as numeric matrix:
            attempting to covert", immediate.=TRUE)
    covar2<-model.matrix(~.,data=as.data.frame(gxe_covar), na.action=na.pass)[,-1]
    print(head(covar2))
    user_input <- readline("Do you want to run using this matrix for covariates (y/n) ")
    if(user_input != 'y') stop('Exiting since you did not press y')
    gxe_covar<-as.matrix(covar2)

  }

  # data entry check
  stopifnot("y must be same length as x" = length(y) == length(x))
  stopifnot("g must be same length as x" = length(g) == length(x))
  stopifnot("covar must have same number of rows  as x" = is.null(covar) ||nrow(covar) == length(x))
  stopifnot("gxe_covar must have same number of rows  as x" = is.null(gxe_covar) ||nrow(gxe_covar) == length(x))
  stopifnot("gxe_interaction must have same number of rows as x" = is.null(gxe_interaction) ||nrow(gxe_interaction) == length(x))
  stopifnot("x_strata must be same length as x" = is.null(x_strata) || length(x_strata) == length(x))
  stopifnot("x_corrected must be same length as x" = is.null(x_corrected) || length(x_corrected) == length(x))

  # coxph checks
  stopifnot(
    "y must be a Surv object with family coxph" = !(family=="coxph" &
                                                      !is.Surv(y))
  )
  stopifnot(
    "cannot use controlsonly option with family coxph" = !(family=="coxph" &
                                                             controlsonly==TRUE)
  )

  ##################### missing data
  #determine who has all data
  complete <- !is.na(x) & !is.na(g) & !is.na(y)
  if (!is.null(covar))           complete <- complete & stats::complete.cases(covar)
  if (!is.null(gxe_covar))       complete <- complete & stats::complete.cases(gxe_covar)
  if (!is.null(gxe_interaction)) complete <- complete & stats::complete.cases(gxe_interaction)
  if (!is.null(x_corrected))     complete <- complete & !is.na(x_corrected)
  if (!is.null(x_strata))        complete <- complete & !is.na(x_strata)

  # drop none complete data and report
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message(n_dropped, " individual(s) removed due to missing data.")
    y <- y[complete]
    x <- x[complete]
    g <- g[complete]
    if (!is.null(covar))           covar           <- covar[complete, , drop = FALSE]
    if (!is.null(gxe_covar))       gxe_covar       <- gxe_covar[complete, , drop = FALSE]
    if (!is.null(gxe_interaction)) gxe_interaction <- gxe_interaction[complete, , drop = FALSE]
    if (!is.null(x_corrected))     x_corrected     <- x_corrected[complete]
    if (!is.null(x_strata))        x_strata        <- x_strata[complete]
  }


###################################### start of function

  # ranking helping function

  calculate_ranked_strata <- function(exp_var, ins_var, prestrat, q) {
    ins_var <- rank(ins_var, ties.method = "random")
    strata1 <- floor((ins_var - 1) / q / prestrat) + 1
    # check GR statistic
    GR_stats<-getGRvalues(X=exp_var, Zstratum=strata1)
    id <- seq_along(exp_var)
    temp <- data.frame(exp_var = exp_var, strata1 = strata1, id = id)
    temp <- arrange(.data = temp, exp_var)
    temp <- group_by(.data = temp, strata1)
    temp <- mutate(.data = temp,
                   x0q = ceiling(rank(exp_var, ties.method = "random") / prestrat))
    temp <- arrange(.data = temp, id)
    list(x0q = temp$x0q, strata = strata1, GR_stats=GR_stats)
  }
  #create strata
 if(is.null(x_strata)){

   # residual approach - calculate the iv-free association
   if (strata_method=="residual"){
   ivf <- iv_free(
       y = y, x = x, g = g,
       covar = covar, gxe_covar=NULL,
       gxe_interaction=NULL,
       q = q, controlsonly=controlsonly
     )
     x0q <- ivf$x0q
   }else if(strata_method=="ranked") {
  # ranked method, rank using x values, unless a residual value is supplied from
     #outside model
    exp_var <- if (!is.null(x_corrected)) x_corrected else x

    ranked<- calculate_ranked_strata(exp_var=exp_var,
                                     ins_var=g,
                                     prestrat = prestrat,
                                     q=q)
    x0q <- ranked$x0q
    GR_stats <- ranked$GR_stats
 }else if(strata_method=="interaction") {
   # ranked method, rank using residual values X - β_(G×E) * (G×E) - iv_free does this if gxe_interaction is not null
   ivf <- iv_free(
     y = y, x = x, g = g,
     covar = NULL, gxe_covar=gxe_covar,
     gxe_interaction=gxe_interaction,
     q = q, controlsonly=controlsonly
   )
   x0 <- ivf$x0
   ranked<- calculate_ranked_strata(exp_var=x0,
                                    ins_var=g,
                                    prestrat = prestrat,
                                    q=q)
   x0q <- ranked$x0q
   GR_stats <- ranked$GR_stats
 }
   quant <- q

 } else {
  if(anyNA(x_strata))(stop("x_strata used; rows missing strata assignment"))
  # relabel to 1..K in case user passed arbitrary labels
  x0q <- match(x_strata, sort(unique(x_strata)))
  quant<- max(x0q)
  }



  # create vector lists
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

  # functions to create values within strata labeled by idx

  # F1: choose boundary probs for xmin/xmax
  xmin_prob <- function(j, quant, strata_bound) if (j == 1) strata_bound[1] else strata_bound[2]
  xmax_prob <- function(j, quant, strata_bound) if (j == quant) strata_bound[3] else strata_bound[4]


  # F2: quantile/mean within stratum
  q_stratum <- function(v, idx, p) quantile(v[idx], probs = p, na.rm = TRUE, names = FALSE, type = 7)
  m_stratum <- function(v, idx) mean(v[idx], na.rm = TRUE)

  # F3: fit outcome-on-g model depending on family and if covariates
  fit_y_on_g <- function(family, y, g, covar, idx) {
    if (!any(idx)) stop("No observations in this stratum.")
    # make dataframe from genetic predictors
    df <- data.frame(g = g)
    # Add covariates (matrix, skip if null)
    if (!is.null(covar)) {
      covar_df <- if (is.vector(covar) && !is.list(covar)) {
        data.frame(V1 = covar)
      } else {
        as.data.frame(covar)
      }
      # Avoid name clashes
      names(covar_df) <- make.names(names(covar_df), unique = TRUE)
      df <- cbind(df, covar_df)
    }
    # Build formula
    if (family == "coxph") {
      # predictors are the columns currently in df (g + covars if present)
      fml <- as.formula(paste("y ~", paste(names(df), collapse = " + ")))
      environment(fml) <- environment()  # ensures y (Surv) is found


    } else {
      df$y <- y
      fml <- if (is.null(covar)) {y ~ g} else {y ~ .}
    }

    # fit relevant model
    model<- switch(
      family,
      binomial = stats::glm(fml, data = df, subset = idx, family = "binomial"),
      coxph    = survival::coxph(fml, data = df, subset = idx),
      gaussian = stats::lm(fml, data = df, subset = idx),
      stop("family must be one of: 'gaussian', 'binomial', 'coxph'")
    )
    model
  }

  # F4: fit exposure-on-g model depending on if covariates
  # and if controls only
  fit_x_on_g <- function(x, y, g, covar, controlsonly, idx, family) {
    if (!any(idx)) stop("No observations in this stratum.")
    if (isTRUE(controlsonly) && family %in% c("binomial", "coxph")) {
      if (is.null(y)) stop("y must be provided when controlsonly=TRUE")
      idx <- idx & (y == 0)
      if (!any(idx)) stop("No controls in this stratum (after applying y==0).")
    }
    # make dataframe
    df <- data.frame(g = g)
    # Add covariates (matrix/dataframe with columns)
    if (!is.null(covar)) {
      covar_df <- if (is.vector(covar) && !is.list(covar)) {
        data.frame(covar1 = covar)
      } else {
        as.data.frame(covar)
      }
      # Avoid name clashes
      names(covar_df) <- make.names(names(covar_df), unique = TRUE)
      df <- cbind(df,covar_df)
    }
    df$x <- x
    # define formula
    fml <- if (is.null(covar)) {
      x ~ g
    } else {
      x ~ .
    }

    # fit relevant model
    stats::lm(fml, data = df, subset = idx)
  }

  # F5: extract coefficient + SE for g from model (and stop nicely if missing)
  coef_se_g <- function(model) {
    coef_temp <- coef(summary(model))
    if (inherits(model, "coxph")) { #gives summary as matrix
     if("g" %in% rownames(coef_temp)) {
        b <- unname(coef_temp["g", "coef"])
        se <- unname(coef_temp["g", "se(coef)"])
      } else {
        b<- NA
        se<-NA
        warning("the regression coefficient in one of the quantiles is missing")
      }
      }else{
        if("g" %in% rownames(coef_temp)) {
          b <- unname(coef_temp["g", "Estimate"])
          se <- unname(coef_temp["g", "Std. Error"])
        }else {
          b<- NA
          se<-NA
          warning("the regression coefficient in one of the quantiles is missing")
        }

      }
    return(list(b=b, se=se))
    }


#use the ivfree quantiles

  for (j in seq_len(quant)) {
    # index the subset
    idx <- (x0q == j)

    # Describe x in stratum
    xmin[j]  <- q_stratum(x, idx, xmin_prob(j, quant, strata_bound))
    xmax[j]  <- q_stratum(x, idx, xmax_prob(j, quant, strata_bound))
    xmean[j] <- m_stratum(x, idx)

    # Fit exposure model (possibly in controls-only)
    mod_x <- fit_x_on_g(x=x, g=g,
                        covar=covar,
                        idx= idx, y = y,
                        controlsonly = controlsonly,
                        family = family)
    x_g <- coef_se_g(mod_x)
    bx[j]   <- x_g$b
    bxse[j] <- x_g$se

    # Fit outcome model
    mod_y <- fit_y_on_g(family =family,
                        y = y,
                        g=g,
                        covar=covar, idx = idx)
    y_g <- coef_se_g(mod_y)
    by[j]   <- y_g$b
    byse[j] <- y_g$se

    if (extra_statistics) {
      strata_stats[[j]] <- list(
        strata = j,
        xmin = q_stratum(x, idx, 0),
        xmax = q_stratum(x, idx, 1),
        xmean = xmean[j],
        ymin = if (family == "coxph") NA_real_ else q_stratum(y, idx, 0),
        ymax = if (family == "coxph") NA_real_ else q_stratum(y, idx, 1),
        x_fstat = (bx[j] / bxse[j])^2
      )
    }

    mod_x <- NULL
    mod_y <- NULL
  }
  # output data

  output <- data.frame(bx, by, bxse, byse, xmean, xmin, xmax)
  names(output) <- c("bx", "by", "bxse", "byse", "xmean", "xmin", "xmax")

  final_output_list= list(summary=output)
  if (extra_statistics) {
    stats<- dplyr::bind_rows(strata_stats)
    final_output_list[["strata_statistics"]]<- stats
  }
  if (strata_method%in%c("ranked", "interaction")){
    if(is.null(x_strata)){
    final_output_list[["GR_max"]]<- GR_stats[1]
    }

    if (report_GR==TRUE && is.null(x_strata)){
      final_output_list[["GR_results"]]<-GR_stats
    }


    if (report_het==TRUE){
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


    p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
    final_output_list[["Heterogeneity_results"]]<- p_heterogeneity
    }
  }
    if(return_assignment){
      final_output_list[["Strata_assignment"]] <- x0q
    }
    final_output_list[["Dropped"]]<-n_dropped


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
#' X, the exposure; and a variety of Y, the outcome values, including
#' \code{y.bin}, a binary outcome generated from a logistic model of X.
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
  lp <- beta1 * data$X + confound * data$u
  data$y.bin <- rbinom(N, 1, stats::plogis(lp - mean(lp)))
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

