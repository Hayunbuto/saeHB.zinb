
# saeHB.zinb

<!-- badges: start -->
<!-- badges: end -->

We designed this package to provide function and datasets for area level of Small Area Estimation under Zero Inflated Negative Binomial Model using Hierarchical Bayesian (HB) Method. This package provides model using Univariate Zero Inflated Negative Binomial Distribution for variable of interest. The 'rjags' package is employed to obtain parameter estimates. Model-based estimators involves the HB estimators which include the mean, the variation of mean, and the quantile of mean. For the reference, see Rao, J.N.K & Molina (2015).

## Author
Hayun , Azka Ubaidillah

## Maintaner
Hayun <221810327@stis.ac.id>

## Installation

You can install the development version of `saeHB.zinb` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hayunbuto/saeHB.zinb")
```

## Function

* `ZinbHB()` This function gives small area estimator under Zero Inflated Negative Binomial Model and is implemented to variable of interest (y) that assumed to be a Zero Inflated Negative Binomial Distribution. The range of data is (y >= 0)


## Example

This is a basic example of using `ZinbHB()` function to make an estimate based on synthetic data in this package

```{r example}
library(saeHB.zinb)
## For data without any non-sampled area
data(dataZINB)       # Load dataset

## For data with non-sampled area use dataHNBNs
## Fitting model
result <- ZinbHB(y ~ x1 + x2, data=dataZINB)
```

Small Area mean Estimates

```r
result$Est
```

Estimated model coefficient

```r
result$coefficient
```

Estimated random effect variances

```r
result$refVar
```


## References
* Desjardins, C. D. (2013). Evaluating the performance of two competing models of school suspension under simulation the zero-inflated negative binomial and the negative binomial hurdle (thesis). Minnesota (US): Minnesota University. <purl.umn.edu/152995>
* Emille E. O. Ishida, Joseph M. Hilbe, and Rafael S. de Souza (2017). Bayesian Models for Astrophysical Data: Using R, JAGS, Python, and Stan. Cambridge : Cambridge University Press. <bayesianmodelsforastrophysicaldata.com>
* Garray, A. M., Hashimoto, E. M., Ortega, E. M. M., dan Lachos, V. H. (2011). On Estimation and Influence Diagnostics For Zero Inflated Negative Binomial Regression Models. Computational Statistics and Data Analysis, 55 (3), p.1304-1318. <doi.org/10.1016/j.csda.2010.09.019>
* Hilbe, J. M. (2011). Negative Binomial Regression 2nd Edition. New York : Cambridge University Press. <doi.org/10.1017/CBO9780511973420>
* Nadhiroh, I. M. (2009). Zero-Inflated Negative Binomial Models in Small Area Estimation. Bogor: Bogor Agricultural University.
* Ntzoufras, I. (2009). Bayesian Modelling Using WinBUGS. New Jersey :  John Wiley & Sons, Inc. <doi.org/10.1002/9780470434567>
* Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New Jersey : John Wiley & Sons, Inc. <doi.org/10.1002/9781118735855>
* S. Krieg, H. J. Boonstra, and M. Smeets. (2016). Small-area estimation with zero-inflated data – a simulation study. J. Off. Stat., vol. 32, no. 4, pp. 963–986, 2016. <doi.org/10.1515/jos-2016-0051>
# saeHB.zinb
