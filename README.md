[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/midasml)](https://cran.r-project.org/package=midasml)
[comment]: [![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/midasml)](https://cran.rstudio.com/web/packages/midasml/index.html) 
[comment]: [![Downloads](http://cranlogs.r-pkg.org/badges/midasml)](http://www.r-pkg.org/pkg/midasml)

# midasML

midasML - estimation and prediction for high-dimensional mixed frequency time series data.

## About

The midasML package implements estimation and prediction methods for high dimensional time series regression models under mixed data sampling data structures using structured-sparsity penalties and orthogonal polynomials. For more information on the midasML approach see [1]. The package also allows to estimate and predict using single-variate MIDAS regressions. Note that such regressions are also implemented in ```midasr``` package. Functions implemented in this package allows to directly compare low-dimensional and high-dimensional MIDAS regression models.

The core of the midasML method is the sparse-group LASSO (sg-LASSO) estimator proposed by [2], and studied for high-dimensional time series data by [1, 3]. The sg-LASSO consists of group structures that are present in high-dimensional ARDL-MIDAS model, hence it is a natural estimator for such model. 

The main algorithm for solving sg-LASSO estimator is taken from [2]. 

Functions that compute MIDAS data structures were inspired by MIDAS Matlab toolbox (v2.3) written by Eric Ghysels, see [4].

## Whats new

ARDL-MIDAS model with Beta density and exponential Almon specifications use analytical gradients derived in [5]. Also, MIDAS Logit model is available for linear in parameters specifications (soon: exponential Almon with analytical gradients)

## Main functions

### Estimation and prediction functions
  - ```midasml_forecast``` - midasML estimation and prediction function.
  - ```midas_ardl``` - ARDL-MIDAS single-variate estimation and prediction function (accomodates different weight functions and loss functions, e.g. quantile regression loss).
  - ```midas_dl``` - DL-MIDAS single-variate estimation and prediction function (accomodates different weight functions and loss functions, e.g. quantile regression loss).
### Estimation only functions
  - ```reg_sgl``` - sg-LASSO regression estimation (currently only for mse loss).
  - ```panel_sgl``` - panel sg-LASSO regression estimation (currently only for mse loss).
### Data handling functions
  - ```qtarget.sort_midasml``` - transforms data into format suitable for midasML technique, creating in-sample and out-of-sample observations for quarterly target variable. Output could be directly inputed into ```midasml_forecast``` (note: currently does not handle real-time data vintages. in case real-time experiment is considered for a specific application, this function can help to setup up the data for each quarter prediction separately. future updates will contain functions capable of handling real-time data vintages.)
  - ```mixed_freq_data``` transforms data into MIDAS regression format creating in-sample and out-of-sample observations. Output is subsequenlty used in ```midas_ardl``` & ```midas_dl```
  
## Run to install the package

```{r }
# CRAN version
install.packages("midasml")

# Development version
# install.packages("devtools")
library(devtools)
install_github("jstriaukas/midasml")
```

## References

[1] Babii, A., Ghysels, E., & Striaukas, J. (2020). Machine learning time series regressions with an application to nowcasting. <https://arxiv.org/abs/2005.14057>

[2] Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013). A sparse-group lasso. Journal of computational and graphical statistics, 22(2), 231-245. Related CRAN R package. https://CRAN.R-project.org/package=SGL 

[3] Babii, A., Ghysels, E., & Striaukas, J. (2020). Inference for high-dimensional regressions with heteroskedasticity and autocorrelation. <https://arxiv.org/abs/1912.06307>.

[4] Ghysels, E. et. al. Mathworks Matlab toolbox. https://www.mathworks.com/matlabcentral/fileexchange/45150-midas-matlab-toolbox

[5] Kostrov (2020) Estimating MIDAS regressions via MIDAS-NLS with revised optimization. Working paper. 
