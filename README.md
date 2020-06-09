# midasML

midasML - estimation and prediction for high-dimensional mixed frequency time series data.

## About
The midasML package constains routines for estimation and prediction methods for high dimensional time series regression models under mixed data sampling data structures using structured-sparsity penalties and orthogonal polynomials. For more information on the midasML approach see [1]. 

The core of the method is the sparse-group LASSO (sg-LASSO) estimator proposed by [2], and studied for high-dimensional time series data by [1]. The main algorithm for solving sg-LASSO estimator is taken from [2]. 

Functions that compute MIDAS data structures were inspired by MIDAS Matlab toolbox (v2.3) written by Eric Ghysels, see [3].

## Main functions












## References

[1] Babii, A., Ghysels, E., & Striaukas, J. (2020). Machine learning time series regressions with an application to nowcasting. <https://arxiv.org/abs/2005.14057>

[2] Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013). A sparse-group lasso. Journal of computational and graphical statistics, 22(2), 231-245. Related CRAN R package. https://CRAN.R-project.org/package=SGL 

[3] Ghysels, E. et. al. Mathworks Matlab toolbox. https://www.mathworks.com/matlabcentral/fileexchange/45150-midas-matlab-toolbox
