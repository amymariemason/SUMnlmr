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
#'   \item{xmean}{Average ldl in stratum}
#'   \item{xmin}{Minimum ldl in stratum}
#'   \item{xmax}{Maximum ldl in stratum}
#'   ...
#' }
"LDL_CAD_covar"

#' Artificial genetic & phenotype data, for purposes of package tests
#'
#' @format A data frame with 10000 rows and 1 variables:
#' \describe{
#'  \item{g}{gene count of 0, 1 or 2. Distibuted binomial with n=2 and p =0.3}
#'  \item{u}{"unmeasured" confounder N(0,1)}
#'   \item{errorX}{Error term for X,  ~exp(1)}
#'   \item{errorY}{Error term for Y, ~N(0.1)}
#'   \item{X}{Exposure. X= 2+ 0.25*g + u + errorX}
#'   \item{linear.Y}{Linear outcome. Y =  X + 0.8u + errorY}
#'   \item{quadratic.Y}{Quadratic outcome. Y= 2X^2 X + 0.8u + error Y}
#'   \item{sqrt.Y}{Square root outcome. \eqn{Y = \sqrt{X} + 0.8*u+ errorY}}
#'   \item{log.Y}{Log outcome. \eqn{Y = \log(X) + 0.8*U +errorY}}
#'   \item{threshold.Y}{\eqn{X + 0.8* U +errorY} if \eqn{X>2}
#'   and \eqn{0.8U + errorY} otherwise}
#'   \item{fracpoly.Y}{Fractional polynomial \eqn{Y = X + 2/X + 0.8U + errorY}}
#' }
"generated_data"

