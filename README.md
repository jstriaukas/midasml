# midasML

midasML - estimation and prediction for high-dimensional mixed frequency time series data.

## About
The midasML package implements estimation and prediction methods for of high dimensional time series regression models under mixed data sampling data structures using structured-sparsity penalties and orthogonal polynomials. For more information on the midasML approach see [1]. Note that single-variate MIDAS regressions are also implemented in ```midasr``` package. Functions for single-variate MIDAS regressions implmented in this package allows to directly compare low-dimensional and high-dimensional MIDAS regression models.

The core of the method is the sparse-group LASSO (sg-LASSO) estimator proposed by [2], and studied for high-dimensional time series data by [1,3]. The main algorithm for solving sg-LASSO estimator is taken from [2]. 

Functions that compute MIDAS data structures were inspired by MIDAS Matlab toolbox (v2.3) written by Eric Ghysels, see [4].

## Main functions

### Estimation and prediction functions
  - ```midasml_forecast``` - midasML estimation and prediction function.
  - ```midas_ardl``` - ARDL-MIDAS single-variate estimation and prediction function (accomodates different weight functions and loss function, e.g. quantile regression loss).
  - ```midas_dl``` - DL-MIDAS single-variate estimation and prediction function (accomodates different weight functions and loss function, e.g. quantile regression loss).
### Estimation only functions
  - ```reg_sgl``` - sg-LASSO regression estimation.
  - ```panel_sgl``` - panel sg-LASSO regression estimation.
### Data handling functions
  - ```qtarget.sort_midasml``` - transforms data into format suitable for midasML technique, creating in-sample and out-of-sample observations. Output could be directly inputed into ```midasml_forecast``` (note: currently does not handle real-time data vintages. in case it real-time experiment is needed, this function could help to setup up data for each quarter prediction separately. future updates will contain functions capable of handling real-time data vintages.)
  - ```mixed_freq_data``` transforms data into MIDAS regression format creating in-sample and out-of-sample observations. Output is subsequenlty used in ```midas_ardl``` & ```midas_dl```
  

## References

[1] Babii, A., Ghysels, E., & Striaukas, J. (2020). Machine learning time series regressions with an application to nowcasting. <https://arxiv.org/abs/2005.14057>

[2] Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013). A sparse-group lasso. Journal of computational and graphical statistics, 22(2), 231-245. Related CRAN R package. https://CRAN.R-project.org/package=SGL 

[3] Babii, A., Ghysels, E., & Striaukas, J. (2020). Inference for high-dimensional regressions with heteroskedasticity and autocorrelation. <https://arxiv.org/abs/1912.06307>.

[4] Ghysels, E. et. al. Mathworks Matlab toolbox. https://www.mathworks.com/matlabcentral/fileexchange/45150-midas-matlab-toolbox
