
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SUMnlmr

<!-- badges: start -->

<!-- badges: end -->

The goal of SUMnlmr is to allow investigations of potentially non-linear
relationships between an exposure and an outcome via a Mendelian
randomization framework, without requiring full access to individual
level genetic data.

It is based on the existing package for individual data by James Staley:
nlmr (available from <https://github.com/jrs95/nlmr> ).

The core concept is to split the process into two distinct halfs: one
requiring individual level data, which is converted into a
semi-summarized form (create\_nlmr\_summary) by dividing the population
into strata based on the IV-free exposure. Associations with the
exposure and the outcome are estimated in each stratum. In the second
half, this semi-summarized form can then be shared, without compromising
patient privacy, and investigated seperately using two IV methods: a
fractional polynomial method (frac\_poly\_summ\_mr) and a piecewise
linear method (piecewise\_summ\_mr). Both methods calculate a localised
causal effect (LACE). The piecewise method fits a continuous piecewise
linear function to these estimates, while the fractional polynomial
method fits the best 1 or 2 term fractional polynomial.

## Functions

*create\_nlmr\_summary* - prepares individual level data into
semi-summarised form, ready to fit nlmr models. *fracpoly\_summ\_mr* -
this method performs IV analysis using fractional polynomials
*piecewise\_summ\_mr* - this method performs IV analysis using piecewise
linear function

## Installation

You can install the released version of SUMnlmr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("amymariemason/SUMnlmr")
```

## Example 1: Summarizing data

This is a basic example which shows you how to create the
semi-summarized data form. First we create some practise data:

``` r
library(SUMnlmr)
## create some data to practise on
test_data<-create_ind_data(N=10000, beta2=2, beta1=1)
# this creates quadratic.Y  = x + 2x^2 + errorY 
head(test_data)
#>   g         u    errorX      errorY        X linear.Y quadratic.Y    sqrt.Y
#> 1 1 0.9005680 2.0223800 -0.09563664 5.172948 5.797766    59.31655 2.8992293
#> 2 0 0.6494483 0.8381424 -1.44817082 3.487591 2.558979    26.88556 0.9388971
#> 3 0 0.4838572 1.3767751  1.51179578 3.860632 5.759514    35.56848 3.8637308
#> 4 0 0.5454714 2.7137997  2.10655517 5.259271 7.802203    63.12207 4.8362424
#> 5 1 0.6899548 0.6736310 -0.60041173 3.613586 3.565138    29.68114 1.8524955
#> 6 0 0.6261793 1.1058177  0.45153692 3.731997 4.684477    32.54008 2.8843181
#>      log.Y threshold.Y fracpoly.Y
#> 1 2.268260    5.797766   9.084651
#> 2 0.320599    2.558979   5.057401
#> 3 3.249713    5.759514   8.461176
#> 4 4.202925    7.802203  11.122188
#> 5 1.236253    3.565138   6.134539
#> 6 2.269424    4.684477   7.318364
```

Then we use create\_nlmr\_summary to summarise it.

``` r

## create the summarized form
## 
summ_data<-create_nlmr_summary(y = test_data$quadratic.Y,
                                x = test_data$X,
                                g = test_data$g,
                                covar = NULL,
                                family = "gaussian",
                                q = 10)

head(summ_data)
#> $summary
#>           bx       by        bxse       byse   x0mean    xmean     xmin
#> 1  0.2429131 2.733407 0.005603564 0.07764084 2.462284 2.420238 2.019524
#> 2  0.2432329 3.012317 0.003179961 0.06079895 2.747566 2.725714 2.604142
#> 3  0.2411304 3.148553 0.002410041 0.05817655 2.939542 2.931291 2.836433
#> 4  0.2469906 3.317581 0.002315272 0.05987202 3.111563 3.108443 3.021684
#> 5  0.2474042 3.528543 0.002752154 0.06212443 3.281760 3.291479 3.202956
#> 6  0.2515016 3.808396 0.003223459 0.07287511 3.484876 3.501065 3.388753
#> 7  0.2453415 4.083432 0.004426916 0.08824239 3.745661 3.760102 3.621663
#> 8  0.2460483 4.300287 0.006051686 0.11896137 4.111155 4.121576 3.913532
#> 9  0.2387040 4.776967 0.009525743 0.19240588 4.637730 4.649872 4.340167
#> 10 0.2381705 5.982458 0.049521500 1.41185388 6.031824 6.044180 5.041532
#>         xmax
#> 1   2.603140
#> 2   2.836406
#> 3   3.021423
#> 4   3.202905
#> 5   3.388672
#> 6   3.621336
#> 7   3.913485
#> 8   4.339599
#> 9   5.039660
#> 10 11.770173
```

If we have co-variants we want to adjust for in our analysis, we need to
include them at this stage.

``` r

## create the summarized form
summ_covar<-create_nlmr_summary(y = test_data$quadratic.Y,
                                x = test_data$X,
                                g = test_data$g,
                                covar = matrix(data=c(test_data$linear.Y,
                                                      test_data$sqrt.Y),ncol=2),
                                family = "gaussian",
                                q = 10)

head(summ_covar)
#> $summary
#>             bx         by         bxse        byse   x0mean    xmean     xmin
#> 1  0.017925800 -1.8303089 1.896658e-03 0.274480400 3.636208 2.420238 2.019524
#> 2  0.007667441 -0.5475614 3.009037e-04 0.022075236 2.993669 2.725714 2.604142
#> 3  0.008125296 -0.6108139 1.656974e-04 0.012628241 3.141104 2.931291 2.836433
#> 4  0.007688979 -0.5970659 1.287311e-04 0.010121586 3.294468 3.108443 3.021684
#> 5  0.007970408 -0.6390177 1.111764e-04 0.009343073 3.476757 3.291479 3.202956
#> 6  0.007841132 -0.6618080 8.662089e-05 0.007803056 3.695819 3.501065 3.388753
#> 7  0.007990362 -0.6958934 8.204351e-05 0.007467348 3.919490 3.760102 3.621663
#> 8  0.007530148 -0.6598165 8.579979e-05 0.007743820 4.027460 4.121576 3.913532
#> 9  0.006924381 -0.6186743 1.086855e-04 0.009737527 4.130589 4.649872 4.340167
#> 10 0.002708821 -0.2437146 2.250764e-04 0.020336782 4.238398 6.044180 5.041532
#>         xmax
#> 1   2.603140
#> 2   2.836406
#> 3   3.021423
#> 4   3.202905
#> 5   3.388672
#> 6   3.621336
#> 7   3.913485
#> 8   4.339599
#> 9   5.039660
#> 10 11.770173
```

These have used a single genetic variant count, but the method works
identically with an genetic score function for g instead.

Once your data is in this format, the output data frame is all you need
to share to fit the fractional polynomial or piecewise linear models
onto the data.

## Example 2: Fitting a fractional polynomial model

Your data needs to be in the semi-summarised form as shown above. We can
then fit a fractional polynomial model:

``` r


model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                  by=by, 
                  bxse=bxse, 
                  byse=byse, 
                  xmean=xmean,
                  family="gaussian",
                  fig=TRUE)
)


summary(model)
#> Call: frac_poly_mr
#> 
#> Number of individuals: NA; Quantiles: 10; 95%CI: Model based SEs
#> 
#> Powers: 2
#> 
#> Coefficients:
#>   Estimate Std. Error 95%CI Lower 95%CI Upper   p.value    
#> 2 2.197367   0.015646    2.166702       2.228 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Fractional polynomial degree p-value: 0.11
#> Fractional polynomial non-linearity p-value: 0
#> Quadratic p-value: 1.21e-75
#> Cochran Q p-value: 0
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0.446
#> Trend p-value: 0.228
```

<img src="man/figures/README-example4-1.png" width="100%" /> This also
produces a graph of the fit with 95% confidence intervals. This is a
ggplot object and can be adjusted with ggplot commands

``` r
library(ggplot2)
f <- function(x) (x + 2*x^2 - mean(summ_data$summary$xmean) -
                    2*mean(summ_data$summary$xmean)^2 )

plot1 <- model$figure+ 
  stat_function(fun = f, colour = "green") +
  ggtitle("fractional polynomial fit from semi-summarized data")

plot1
```

<img src="man/figures/README-example 5-1.png" width="100%" /> There is
also p-values provided in p\_test and p\_het. This is identical to the
testing provided by the nlmr package: \* fp\_d1\_d2 : test between the
fractional polynomial degrees \* fp : fractional polynomial
non-linearity test \* quad: quadratic test \* Q : Cochran Q test and \*
Q: Cochran Q heterogeneity test \* trend: trend test

``` r
model$p_tests
#>       fp_d1_d2 fp         quad Q
#> [1,] 0.1104167  0 1.211889e-75 0

model$p_heterogeneity
#>              Q     trend
#> [1,] 0.4458058 0.2276787
```

## Example 3: Piecewise linear model

We can instead fit a piecewise linear model to the same summarised data

``` r

model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, x0mean, xmean, xmin,xmax, 
                  ci="bootstrap_se",
                  nboot=1000, 
                  fig=TRUE,
                  family="gaussian",
                  ci_fig="ribbon")
)

summary(model2)
#> Call: piecewise_mr
#> 
#> Number of individuals: ; Quantiles: 10; Number of bootstrap replications: 1000
#> 
#> LACE:
#>    Estimate Std. Error 95%CI Lower 95%CI Upper   p.value    
#> 1  11.13532    0.31629    10.49014      11.780 < 2.2e-16 ***
#> 2  12.27154    0.24768    11.80163      12.742 < 2.2e-16 ***
#> 3  12.82654    0.23700    12.36498      13.288 < 2.2e-16 ***
#> 4  13.51512    0.24391    13.01722      14.013 < 2.2e-16 ***
#> 5  14.37454    0.25308    13.85346      14.896 < 2.2e-16 ***
#> 6  15.51460    0.29688    14.94285      16.086 < 2.2e-16 ***
#> 7  16.63504    0.35948    15.92894      17.341 < 2.2e-16 ***
#> 8  17.51846    0.48462    16.55800      18.479 < 2.2e-16 ***
#> 9  19.46035    0.78382    17.92542      20.995 < 2.2e-16 ***
#> 10 24.37127    5.75159    13.03990      35.703 2.492e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Quadratic p-value: 1.32e-75
#> Cochran Q p-value: 0
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0.446
#> Trend p-value: 0.241
```

<img src="man/figures/README-pl example-1.png" width="100%" /> Again the
figure is a ggplot object and can be adjusted similarly.

``` r
plot2 <- model2$figure+ 
  stat_function(fun = f, colour = "green") +
  ggtitle("piecewise linear fit from semi-summarized data")

plot2
```

<img src="man/figures/README-pl2-1.png" width="100%" />

## Example 4: Binary outcome

The functions above can also fit binary outcome data, via a generalised
linear model.

``` r

test_data$y.bin<-stats::rbinom(size=1, p=0.5, n=10000)

# create semo-summ data
summ_bin<-create_nlmr_summary(y = test_data$y.bin,
                                x = test_data$X,
                                g = test_data$g,
                                covar = NULL,
                                family = "binomial",
                                q = 10)

# fit fractional poly model


model3<- with(summ_bin$summary,frac_poly_summ_mr(bx=bx,
                  by=by, 
                  bxse=bxse, 
                  byse=byse, 
                  xmean=xmean,
                  family="gaussian",
                  fig=TRUE)
)

summary(model3)
#> Call: frac_poly_mr
#> 
#> Number of individuals: NA; Quantiles: 10; 95%CI: Model based SEs
#> 
#> Powers: -2
#> 
#> Coefficients:
#>    Estimate Std. Error 95%CI Lower 95%CI Upper p.value
#> -2  -1.0506     1.8704     -4.7165      2.6153  0.5743
#> 
#> Non-linearity tests
#> Fractional polynomial degree p-value: 0.943
#> Fractional polynomial non-linearity p-value: 0.676
#> Quadratic p-value: 0.97
#> Cochran Q p-value: 0.751
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0.571
#> Trend p-value: 0.251
```

<img src="man/figures/README-bin-1.png" width="100%" /> Not
unsurprisingly, we find no evidence of an effect, causal or otherwise,
as the binary outcome was randomly distributed.

If we look instead at the semi-summarised UK Biobank datasets on
LDL-cholesterol and CAD, one with and one without covariates. Here we
can see a potentially non-linear trend in the univariate data, which
becomes a clear linear trend once covariates are included.

``` r

# fit piecewise linear model
model4 <-with(LDL_CAD, piecewise_summ_mr(by, bx, byse, bxse, x0mean, xmean, xmin,xmax, 
                  ci="bootstrap_se",
                  nboot=1000, 
                  fig=TRUE,
                  family="gaussian",
                  ci_fig="ribbon")
)


summary(model4)
#> Call: piecewise_mr
#> 
#> Number of individuals: ; Quantiles: 10; Number of bootstrap replications: 1000
#> 
#> LACE:
#>     Estimate Std. Error 95%CI Lower 95%CI Upper   p.value    
#> 1   0.406261   0.059553    0.291876      0.5206 3.371e-12 ***
#> 2   0.362442   0.073477    0.215760      0.5091 1.279e-06 ***
#> 3   0.311077   0.087016    0.139440      0.4827 0.0003818 ***
#> 4   0.286820   0.096907    0.094762      0.4789 0.0034217 ** 
#> 5   0.152663   0.101981   -0.046355      0.3517 0.1327153    
#> 6   0.141322   0.105426   -0.072826      0.3555 0.1958532    
#> 7   0.128421   0.102452   -0.070173      0.3270 0.2050004    
#> 8   0.191138   0.104945   -0.011314      0.3936 0.0642469 .  
#> 9   0.223650   0.104324    0.015068      0.4322 0.0355888 *  
#> 10  0.268987   0.098653    0.074907      0.4631 0.0065979 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Quadratic p-value: 0.0134
#> Cochran Q p-value: 0.193
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0
#> Trend p-value: 4.47e-19
```

<img src="man/figures/README-bi2-1.png" width="100%" />

``` r


# fit piecewise linear model
model5 <-with(LDL_CAD_covar,piecewise_summ_mr(by, bx, byse, bxse, x0mean, xmean, xmin,xmax, 
                  ci="bootstrap_se",
                  nboot=1000, 
                  fig=TRUE,
                  family="gaussian",
                  ci_fig="ribbon")
)

summary(model5)
#> Call: piecewise_mr
#> 
#> Number of individuals: ; Quantiles: 10; Number of bootstrap replications: 1000
#> 
#> LACE:
#>     Estimate Std. Error 95%CI Lower 95%CI Upper   p.value    
#> 1  0.3076383  0.0628626   0.1851991      0.4301 8.451e-07 ***
#> 2  0.2922888  0.0780399   0.1426731      0.4419 0.0001286 ***
#> 3  0.3579244  0.0915492   0.1813115      0.5345 7.123e-05 ***
#> 4  0.2152398  0.1014575   0.0062028      0.4243 0.0435747 *  
#> 5  0.2541277  0.1062581   0.0500576      0.4582 0.0146555 *  
#> 6  0.4308538  0.1079399   0.2211773      0.6405 5.637e-05 ***
#> 7  0.2584994  0.1068753   0.0468446      0.4702 0.0166749 *  
#> 8  0.2921218  0.1075173   0.0798857      0.5044 0.0069811 ** 
#> 9  0.3296967  0.1062308   0.1228140      0.5366 0.0017869 ** 
#> 10 0.3950466  0.1005297   0.1959352      0.5942 0.0001008 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Quadratic p-value: 0.555
#> Cochran Q p-value: 0.928
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0
#> Trend p-value: 7.45e-19
```

<img src="man/figures/README-bi2-2.png" width="100%" />
