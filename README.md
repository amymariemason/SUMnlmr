
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
#>   g         u    errorX     errorY        X linear.Y quadratic.Y    sqrt.Y
#> 1 0 0.2489799 0.9543541  1.1770164 3.203334 4.579534    25.10223 3.1659864
#> 2 2 0.1755614 2.1280490  0.5508806 4.803610 5.494940    51.64429 2.8830437
#> 3 0 0.1271525 0.1282090 -0.3941085 2.255362 1.962975    12.13629 1.2093996
#> 4 1 0.8511619 0.3516511 -1.7530306 3.452813 2.380712    26.22455 0.7860736
#> 5 0 0.7960857 0.1932762 -1.8675103 2.989362 1.758720    19.63129 0.4983354
#> 6 1 0.7606432 2.7672322 -0.1747528 5.777875 6.211637    72.97932 2.8374829
#>        log.Y threshold.Y fracpoly.Y
#> 1  2.5403925    4.579534   6.907919
#> 2  2.2606975    5.494940   8.633676
#> 3  0.5209238    1.962975   3.589596
#> 4  0.1670882    2.380712   4.859090
#> 5 -0.1355818    1.758720   3.948840
#> 6  2.1877977    6.211637   9.719709
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
                                controlsonly = FALSE,
                                q = 10)

head(summ_data)
#> $summary
#>           bx       by        bxse       byse    xmean     xmin     xmax
#> 1  0.2588236 2.868310 0.005621032 0.08044786 2.475200 2.306240 2.734920
#> 2  0.2639698 3.262953 0.003193135 0.06357513 2.753109 2.530231 2.967550
#> 3  0.2663603 3.490383 0.002404029 0.05902643 2.947150 2.741072 3.137965
#> 4  0.2632726 3.673075 0.002233041 0.05738331 3.122582 2.916324 3.300493
#> 5  0.2654057 3.832009 0.002526511 0.06608632 3.275774 3.073859 3.463344
#> 6  0.2657900 3.968406 0.003074774 0.06710325 3.480791 3.250359 3.686875
#> 7  0.2662545 4.202888 0.004110042 0.08104480 3.721853 3.474006 3.957029
#> 8  0.2667502 4.636122 0.005614201 0.10936277 4.044472 3.771661 4.333903
#> 9  0.2388100 4.658954 0.009509496 0.18909153 4.559325 4.226539 4.902271
#> 10 0.3452849 8.874508 0.045877868 1.41920866 5.919853 5.058444 6.537480
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
#>             bx         by         bxse        byse    xmean     xmin     xmax
#> 1  0.020112729 -2.5202320 2.628699e-03 0.502447680 3.517292 2.306240 7.179041
#> 2  0.007675647 -0.5513467 2.781174e-04 0.020213357 2.998743 2.555001 2.916064
#> 3  0.007540994 -0.5649603 1.703104e-04 0.012959081 3.184772 2.757848 5.578689
#> 4  0.007438911 -0.5801280 1.234567e-04 0.009951645 3.315705 2.916463 5.179518
#> 5  0.007400747 -0.5863771 9.790026e-05 0.008083955 3.476401 3.045057 5.311324
#> 6  0.007495053 -0.6212495 9.193423e-05 0.007956212 3.648358 3.172838 5.156910
#> 7  0.007345388 -0.6334340 7.316545e-05 0.006657355 3.842452 3.328465 4.996101
#> 8  0.007020022 -0.6158171 7.685997e-05 0.006920936 4.000645 3.459216 4.861599
#> 9  0.006524188 -0.5769882 1.046986e-04 0.009336814 4.111890 3.623456 4.749192
#> 10 0.002731519 -0.2439639 2.130224e-04 0.018990326 4.203852 3.781399 4.493341
```

Note: Because the covariants are included as a matrix, lm cannot detect
factor variables and create automatic dummy variables for them. The
easiest way to include factor variables is to make these dummy variables
by hand instead using the
[model.matrix](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/model.matrix.html)
command. e.g.

``` r

## create a factor
test_data$centre<- as.factor(rbinom(nrow(test_data),4, 0.5))

#turn factor into binary contrasts against first factor
dummies<- model.matrix(~centre,data=test_data)[,2:5]

summ_covar2<-create_nlmr_summary(y = test_data$quadratic.Y,
                                 x = test_data$X,
                                 g = test_data$g,
                                 covar = dummies,
                                 family = "gaussian",
                                 q = 10)

head(summ_covar2)
#> $summary
#>           bx       by        bxse       byse    xmean     xmin     xmax
#> 1  0.2559522 2.818961 0.005664330 0.08046629 2.473610 2.306240 2.735941
#> 2  0.2614943 3.248641 0.003156374 0.06340713 2.754592 2.533839 2.968272
#> 3  0.2654447 3.475747 0.002406486 0.05972361 2.946783 2.741906 3.139787
#> 4  0.2627423 3.679040 0.002199282 0.05728084 3.124169 2.917619 3.309920
#> 5  0.2651361 3.830524 0.002558838 0.06664064 3.274785 3.073161 3.466214
#> 6  0.2664208 3.980504 0.003103842 0.06859975 3.482717 3.250688 3.694950
#> 7  0.2662619 4.187337 0.004146708 0.08214068 3.720149 3.471521 3.958958
#> 8  0.2680148 4.670709 0.005631446 0.11091571 4.045421 3.771661 4.336108
#> 9  0.2422790 4.723140 0.009477271 0.18882697 4.558720 4.222653 4.906754
#> 10 0.3485304 8.934549 0.045833562 1.41926455 5.919162 5.054255 6.537480
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
#> 2 2.188809   0.014287    2.160806      2.2168 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Fractional polynomial degree p-value: 0.00382
#> Fractional polynomial non-linearity p-value: 0
#> Quadratic p-value: 1.02e-72
#> Cochran Q p-value: 0
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0.163
#> Trend p-value: 0.939
```

<img src="man/figures/README-example4-1.png" width="100%" /> This also
produces a graph of the fit with 95% confidence intervals. This is a
ggplot object and can be adjusted with ggplot commands

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.0.5
f <- function(x) (x + 2*x^2 - mean(summ_data$summary$xmean) -
                    2*mean(summ_data$summary$xmean)^2 )

plot1 <- model$figure+ 
  stat_function(fun = f, colour = "green") +
  ggtitle("fractional polynomial fit from semi-summarized data")

plot1
```

<img src="man/figures/README-example5-1.png" width="100%" /> There is
also p-values provided in p\_test and p\_het. This is identical to the
testing provided by the nlmr package: \* fp\_d1\_d2 : test between the
fractional polynomial degrees \* fp : fractional polynomial
non-linearity test \* quad: quadratic test \* Q : Cochran Q test and \*
Q: Cochran Q heterogeneity test \* trend: trend test

``` r
model$p_tests
#>        fp_d1_d2 fp         quad Q
#> [1,] 0.00382117  0 1.016544e-72 0

model$p_heterogeneity
#>             Q     trend
#> [1,] 0.162997 0.9391521
```

## Example 3: Piecewise linear model

We can instead fit a piecewise linear model to the same summarised data

``` r

model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                  ci="bootstrap_se",
                  nboot=1000, 
                  fig=TRUE,
                  family="gaussian",
                  ci_fig="ribbon")
)

summary(model2)
#> $model
#>    q nboot
#> 1 10  1000
#> 
#> $coefficients
#>        beta        se      lci      uci          pval
#> 1  10.84182 0.3040819 10.23736 11.44628 9.694168e-271
#> 2  12.33351 0.2403053 11.85092 12.81611  0.000000e+00
#> 3  13.19317 0.2231118 12.75837 13.62797  0.000000e+00
#> 4  13.88372 0.2169010 13.44159 14.32586  0.000000e+00
#> 5  14.48447 0.2497972 13.96163 15.00731  0.000000e+00
#> 6  15.00003 0.2536411 14.48284 15.51723  0.000000e+00
#> 7  15.88634 0.3063382 15.26497 16.50771  0.000000e+00
#> 8  17.52391 0.4133763 16.70458 18.34323  0.000000e+00
#> 9  17.61021 0.7147401 16.13224 19.08818 1.263645e-120
#> 10 33.54443 5.3644143 23.11856 43.97029  2.860621e-10
#> 
#> $p_tests
#>           quad Q
#> 1 1.016544e-72 0
#> 
#> $p_heterogeneity
#>          Q     trend
#> 1 0.162997 0.8732123
#> 
#> $figure
```

<img src="man/figures/README-plexample-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "summary.piecewise_mr"

Again the figure is a ggplot object and can be adjusted similarly.

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
                  family="binomial",
                  fig=TRUE)
)
#> Warning: Fisher scoring algorithm may have gotten stuck at a local maximum.
#>   Setting tau^2 = 0. Check the profile likelihood plot with profile().

summary(model3)
#> Call: frac_poly_mr
#> 
#> Number of individuals: NA; Quantiles: 10; 95%CI: Model based SEs
#> 
#> Powers: -1
#> 
#> Coefficients:
#>    Estimate Std. Error 95%CI Lower 95%CI Upper p.value  
#> -1  -2.2290     1.0791     -4.3441     -0.1139 0.03887 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Non-linearity tests
#> Fractional polynomial degree p-value: 0.497
#> Fractional polynomial non-linearity p-value: 0.149
#> Quadratic p-value: 0.0683
#> Cochran Q p-value: 0.479
#> 
#> Heterogeneity tests
#> Cochran Q p-value: 0.0147
#> Trend p-value: 0.183
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
model4 <-with(LDL_CAD, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                  ci="bootstrap_se",
                  nboot=1000, 
                  fig=TRUE,
                  family="gaussian",
                  ci_fig="ribbon")
)


summary(model4)
#> $model
#>    q nboot
#> 1 10  1000
#> 
#> $coefficients
#>         beta         se         lci       uci         pval
#> 1  0.4062614 0.05955269  0.29040291 0.5221199 6.295310e-12
#> 2  0.3624422 0.07347687  0.21920787 0.5056765 7.063386e-07
#> 3  0.3110766 0.08701566  0.13598922 0.4861640 4.970907e-04
#> 4  0.2868202 0.09690730  0.09966541 0.4739751 2.666680e-03
#> 5  0.1526631 0.10198090 -0.05107052 0.3563968 1.419193e-01
#> 6  0.1413218 0.10542626 -0.06461497 0.3472586 1.786162e-01
#> 7  0.1284211 0.10245230 -0.07345933 0.3303016 2.124693e-01
#> 8  0.1911379 0.10494505 -0.01139742 0.3936733 6.435630e-02
#> 9  0.2236501 0.10432433  0.01842177 0.4288784 3.268478e-02
#> 10 0.2689869 0.09865336  0.07254420 0.4654297 7.278910e-03
#> 
#> $p_tests
#>        quad         Q
#> 1 0.0129998 0.1926024
#> 
#> $p_heterogeneity
#>   Q        trend
#> 1 0 3.752917e-19
#> 
#> $figure
```

<img src="man/figures/README-bi2-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "summary.piecewise_mr"
    
    
    # fit piecewise linear model
    model5 <-with(LDL_CAD_covar,piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                      ci="bootstrap_se",
                      nboot=1000, 
                      fig=TRUE,
                      family="gaussian",
                      ci_fig="ribbon")
    )
    
    summary(model5)
    #> $model
    #>    q nboot
    #> 1 10  1000
    #> 
    #> $coefficients
    #>         beta         se        lci       uci         pval
    #> 1  0.3076383 0.06286262 0.18442139 0.4308552 9.902095e-07
    #> 2  0.2922888 0.07803988 0.14098715 0.4435904 1.528555e-04
    #> 3  0.3579244 0.09154918 0.18321943 0.5326293 5.931395e-05
    #> 4  0.2152398 0.10145750 0.01946529 0.4110144 3.117145e-02
    #> 5  0.2541277 0.10625811 0.04565396 0.4626015 1.688404e-02
    #> 6  0.4308538 0.10793986 0.22351069 0.6381970 4.644522e-05
    #> 7  0.2584994 0.10687532 0.04024621 0.4767526 2.026381e-02
    #> 8  0.2921218 0.10751727 0.08258085 0.5016627 6.286651e-03
    #> 9  0.3296967 0.10623077 0.12604819 0.5333452 1.507990e-03
    #> 10 0.3950466 0.10052968 0.19464311 0.5954501 1.116999e-04
    #> 
    #> $p_tests
    #>        quad         Q
    #> 1 0.5564354 0.9280444
    #> 
    #> $p_heterogeneity
    #>   Q        trend
    #> 1 0 5.974956e-19
    #> 
    #> $figure

<img src="man/figures/README-bi2-2.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "summary.piecewise_mr"
