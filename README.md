
Overview
--------

psindex is a wrapper package for lavaan that evaluates parameter estimate stability. The psindex allows for a comparison between model maximum likelihood estimates and other alternative estimates that similarly describe the model (i.e., fungible estimates). The ranges of these alternative estimates allow for an examination of uncertainty in ones results.

Installation
------------

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("nietsnel/psindex")
```

Usage
-----

To start using psindex, the following example allows for the examination of uncertainty in the Political Democracy dataset referenced in Bollen's 1989 text (p. 324)

``` r
library(psindex)
model_original <-
'# measurement model
ind60 =~x1 + x2 + x3
dem60 =~y1 + y2 + y3 + y4
dem65 =~y5 + y6 + y7 + y8
# regressions
dem60∼ind60
dem65∼ind60 + dem60
# residual correlations
y1~~y5
y2~~y4 + y6
y3~~y7
y4~~ y8
y6~~y8
'

ps_index(model = model_original, data_set = PoliticalDemocracy,
         RMSEA_pert = .005,
         plot_fpe =  TRUE, output_long = FALSE,
         frac_plot = .3, iterations_bin = 100000,
         control_genSA = list(threshold.stop = 1e-13, max.time = 600))
```

Reference
---------

Prendez, J. Y., & Harring, J. R. (2019). Measuring Parameter Uncertainty by Identifying &nbsp Fungible Estimates in SEM. Structural Equation Modeling: A Multidisciplinary Journal, 0(0), 1–12. <https://doi.org/10.1080/10705511.2019.1608550>
