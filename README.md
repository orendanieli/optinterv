
optinterv: optimal intervention
===============================

[![CRAN status](https://www.r-pkg.org/badges/version/optinterv)](https://cran.r-project.org/package=optinterv) [![metacran downloads](https://cranlogs.r-pkg.org/badges/optinterv)](https://cran.r-project.org/package=optinterv)

An R package that implemnts the method proposed by Danieli, Devi and Fryer (2019), to identify the factors with the greatest potential to increase a pre-specified outcome, using observational data.

Installation
------------

You can install the stable version from CRAN:

``` r
install.packages("optinterv")
```

Or the development version from GitHub:

``` r
#install.packages("devtools")
devtools::install_github("eladg9/optinterv")
```

Useage Example
--------------

``` r
library(optinterv)
#generate data
n <- 1000
p <- 10
features <- matrix(rnorm(n*p), ncol = p)
men <- matrix(rbinom(n, 1, 0.5), nrow = n)
outcome <- 2*(features[,1] > 1) + men*pmax(features[,2], 0) + rnorm(n)
outcome <- as.vector(outcome)
#find the optimal intervention using the non-parametric method:
imp_feat <- optint(Y = outcome, X = features, control = men, 
                   method = "non-parametric", lambda = 10, plot = TRUE)
#by default, only the significant features are displayed 
#(see ?plot.optint for further details).
#for customized variable importance plot, use plot():
plot(imp_feat, plot.vars = 10)
#show summary of the results using summary():
summary(imp_feat)
#we can look on the new features distribution more deeply, using plot_change():
plot_change(imp_feat, plot.vars = "sig")
#we can explore how the optimal intervention varies between genders using optint_by_group():
men <- as.vector(men)
imp_feat_by_gender <- optint_by_group(Y = outcome, X = features,
                                      group = men, 
                                      method = "non-parametric",
                                      lambda = 10)
#by default, only the significant features are displayed 
#(see ?plot.optint_by_group for further details). 
#for customized variable importance plot, use plot():
plot(imp_feat_by_gender, plot.vars = 10)
```

References
----------
