
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
test_data<-create_ind_data(N=10000)
head(test_data)
#>   g          u    errorX     errorY        X  linear.Y quadratic.Y   sqrt.Y
#> 1 1 0.40805622 3.0064403 -0.2995174 5.664496 17.020417   241.62606 7.166988
#> 2 1 0.21534885 0.6446562 -0.8621351 3.110005  8.640159    76.34508 4.600706
#> 3 0 0.74804529 1.3708998  1.4614011 4.118945 14.416673   133.17663 8.148393
#> 4 0 0.72256493 0.3073310 -1.5512226 3.029896  8.116517    72.37840 4.248808
#> 5 0 0.20230427 0.3063129  0.1892006 2.508617  7.876896    51.92902 5.102628
#> 6 0 0.06207741 0.7464886 -0.3607367 2.808566  8.114623    63.33092 4.716558
#>      log.Y threshold.Y fracpoly.Y
#> 1 5.229582  0.02692756   29.15994
#> 2 2.714017 -0.68985598   16.58253
#> 3 6.306629  2.05983738   24.32585
#> 4 2.352414 -0.97317070   15.87622
#> 5 3.110239  0.35104406   14.31502
#> 6 2.786947 -0.31107474   15.34334
```

Then we use create\_nlmr\_summary to summarise it.

``` r

## create the summarized form
summ_data<-create_nlmr_summary(y = test_data$quadratic.Y,
                                x = test_data$X,
                                g = test_data$g,
                                covar = NULL,
                                family = "gaussian",
                                q = 10)

head(summ_data)
#> $summary
#>           bx        by        bxse      byse   x0mean    xmean     xmin
#> 1  0.2404186  9.096816 0.005688286 0.2201619 2.453012 2.408442 2.005664
#> 2  0.2507362 10.477813 0.003156837 0.1415596 2.731423 2.716669 2.593070
#> 3  0.2468633 11.080452 0.002598484 0.1259060 2.943620 2.928386 2.830601
#> 4  0.2489448 11.767321 0.002173975 0.1126754 3.116760 3.111005 3.024371
#> 5  0.2471969 12.361046 0.002461135 0.1300851 3.282905 3.292795 3.202491
#> 6  0.2500613 13.210698 0.003421844 0.1842867 3.481310 3.501868 3.391309
#> 7  0.2423650 13.604444 0.004109931 0.2323565 3.732510 3.745949 3.617463
#> 8  0.2392425 14.521088 0.005562188 0.3379124 4.073458 4.085273 3.895587
#> 9  0.2277081 15.499179 0.009830862 0.6667058 4.593274 4.606786 4.302539
#> 10 0.2511986 22.333190 0.046459969 4.6367561 5.973015 5.984113 4.990727
#>         xmax
#> 1   2.593029
#> 2   2.830485
#> 3   3.024272
#> 4   3.202366
#> 5   3.391268
#> 6   3.617089
#> 7   3.895199
#> 8   4.302488
#> 9   4.988986
#> 10 11.857938
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
#>             bx         by         bxse       byse   x0mean    xmean     xmin
#> 1  0.015396537 -5.1750704 1.958950e-03 1.01136364 3.509843 2.408442 2.005664
#> 2  0.008148862 -2.0676958 2.708160e-04 0.07128303 3.070771 2.716669 2.593070
#> 3  0.007727882 -2.0590011 1.736100e-04 0.04809538 3.230546 2.928386 2.830601
#> 4  0.007832965 -2.1140258 1.262623e-04 0.03613632 3.320888 3.111005 3.024371
#> 5  0.007999900 -2.2355057 9.298827e-05 0.02866328 3.477625 3.292795 3.202491
#> 6  0.007820785 -2.2789625 8.806538e-05 0.02779696 3.687261 3.501868 3.391309
#> 7  0.007568582 -2.2464194 7.691978e-05 0.02418114 3.855189 3.745949 3.617463
#> 8  0.007648772 -2.2924873 8.305179e-05 0.02565449 3.970591 4.085273 3.895587
#> 9  0.006961814 -2.1189077 1.087843e-04 0.03361986 4.093882 4.606786 4.302539
#> 10 0.002390358 -0.7288713 2.206655e-04 0.06719389 4.164691 5.984113 4.990727
#>         xmax
#> 1   2.593029
#> 2   2.830485
#> 3   3.024272
#> 4   3.202366
#> 5   3.391268
#> 6   3.617089
#> 7   3.895199
#> 8   4.302488
#> 9   4.988986
#> 10 11.857938
```

These have used a single genetic variant count, but the method works
identically with an genetic score function for g instead.

Once your data is in this format, the output data frame is all you need
to share to fit the fractional polynomial or piecewise linear models
onto the data.

## Example 2: Fitting a fractional polynomial model

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
