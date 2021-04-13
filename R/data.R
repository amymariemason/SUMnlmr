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
