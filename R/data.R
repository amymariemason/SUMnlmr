#' Hunt dataset containing the allele score-associations with bmi and all cause
#' mortality, by 100 quantiles of residual BMI. As used in the paper
#' "Body mass index and all cause mortality in HUNT and UK Biobank studies:
#' linear and non-linear mendelian randomisation analyses"
#' https://doi.org/10.1136/bmj.l1042
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{beta_bmi}{linear association between BMI and allele score, adjusted
#'   for age, sex and age-squared}
#'   \item{se_bmi}{standard error for the beta_bmi term from the linear
#'   regression}
#'   \item{beta_acm}{Cox proportional hazards regression association between all
#'    cause mortality and allele score, adjusted for age and sex}
#'   \item{se_acm}{standard error for the beta_acm term from the Cox
#'   regression}
#'   \item{mean_bmi}{Average BMI in stratum}
#'   ...
#' }
"bmi_acm"

#' UK Biobank dataset containing the allele score-associations with
#' ldl-cholesterol and CAD diagnosis, by 10 quantiles of residual ldl.
#'
#' @format A data frame with 10 rows and 8 variables:
#' \describe{
#'   \item{bx}{linear association between ldl-cholesterol and allele score,
#'   unadjusted}
#'    \item{by}{logistic regression coefficent between CAD diagnosis and allele
#'    score, unadjusted}
#'   \item{bxse}{standard error for the bx term from the linear
#'   regression}
#'   \item{byse}{standard error for the by term from the logistic
#'   regression}
#'   \item{x0mean}{Average residual ldl in stratum}
#'   \item{xmean}{Average ldl in stratum}
#'   \item{xmin}{Minimum ldl in stratum}
#'   \item{xmax}{Maximum ldl in stratum}
#'   ...
#' }
"LDL_CAD"

#' UK Biobank dataset containing the allele score-associations with
#' ldl-cholesterol and CAD diagnosis adjusted for age, sex and first 10
#' principle components, by 10 quantiles of residual ldl.
#'
#' @format A data frame with 10 rows and 8 variables:
#' \describe{
#'   \item{bx}{linear association between ldl-cholesterol and allele score,
#'   adjusted}
#'    \item{by}{logistic regression coefficent between CAD diagnosis and allele
#'    score, adjusted}
#'   \item{bxse}{standard error for the bx term from the linear
#'   regression}
#'   \item{byse}{standard error for the by term from the logistic
#'   regression}
#'   \item{x0mean}{Average residual ldl in stratum}
#'   \item{xmean}{Average ldl in stratum}
#'   \item{xmin}{Minimum ldl in stratum}
#'   \item{xmax}{Maximum ldl in stratum}
#'   ...
#' }
"LDL_CAD_covar"
