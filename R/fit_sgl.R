#' Panel sg-LASSO regression model
#' 
#' @description 
#' Fits panel data regression model, random and fixed effects, under sg-LASSO penalty function. Options include random effects and fixed effects models, cross-validation and information criteria for \eqn{\lambda} penalty parameter selection. 
#' @details
#' \ifelse{html}{\out{The sequence of panel data models implied by <code>lambdas</code> vector is fit by block coordinate-descent. The objective function is <br><br> <center> RSS(&alpha;,&beta;)/NT + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where RSS(&alpha;,&beta;) is either random or fixed effects model fit. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on coefficients of time-varying covariates &beta; and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions. There are additional options to apply sg-LASSO structures on fixed effects coefficient vector &alpha; (see description of input variables for more details).}}{The sequence of panel data models implied by \code{lambdas} vector is fit by block coordinate-descent. The objective function is \deqn{RSS(\alpha,\beta)/NT + 2\lambda * \Omega_\gamma(\beta),} where \eqn{RSS(\alpha,\beta)} is either random or fixed effects model fit. The penalty function \eqn{\Omega_\gamma(.)} is applied on coefficients of time-varying covariates \eqn{\beta} and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_{2,1},} a convex combination of LASSO and group LASSO penalty functions. There are additional options to apply sg-LASSO structures on fixed effects coefficient vector \eqn{\alpha} (see description of input variables for more details).}     
#' @usage 
#' panel_sgl(X, Z = NULL, y, index, entity_indices, gamma_w = NULL, l1_factor = NULL, 
#'   l21_factor = NULL, dummies_index = NULL, full_est = NULL, 
#'   regress_choice = c("re", "fe"), method_choice = c("ic", "cv", "initial"), 
#'   nlam = 100, lambdas = NULL, min_frac = NULL, nfolds = 10, 
#'   lambda_choice = c("min", "1se"), ic_choice = c("bic", "aic", "aicc"),
#'   num_cores = NULL, verbose = FALSE, thresh = NULL, 
#'   outer_thresh = NULL, inner_iter = NULL, outer_iter = NULL)
#' @param X NT by p data matrix, where n, t and p respectively denote the number of individuals, sample size and the number of regressors.
#' @param Z dummies matrix for random effects or fixed effects panel data model. If left unspecified, it is computed based on \code{regress_choice} choice.
#' @param y NT by 1 vector of outcome.
#' @param index p by 1 vector indicating group membership of each covariate.
#' @param entity_indices  NT by 1 vector of individual indices.
#' @param gamma_w sg-LASSO mixing parameter. \code{gamma_w = 1} is LASSO and \code{gamma_w = 0} group LASSO.
#' @param l1_factor \ifelse{html}{\out{&#8467;<sub>1</sub> norm penalty factor for random or fixed effects (default value is zero which means &alpha; is left unpenalized in &#8467;<sub>1</sub> norm).}}{\eqn{\ell_1} norm penalty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1} norm).}     
#' @param l21_factor \ifelse{html}{\out{&#8467;<sub>2,1</sub> norm penalty factor for random or fixed effects (default value is zero which means &alpha; is left unpenalized in &#8467;<sub>2,1</sub> norm).}}{\eqn{\ell_1} norm penalty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1} norm).}     
#' @param dummies_index vector indicating group membership of \eqn{\alpha} (default - no grouping).
#' @param full_est pre-estimated parameters based on full sample and \code{regress_choice} for a sequence of \eqn{\lambda}'s.
#' @param regress_choice choose between `re` and `fe`. `re` computes random effects regression with sg-LASSO penalty (default). `fe` computes fixed effects regression with sg-LASSO penalty.
#' @param method_choice choose between `initial`, `ic` and `cv`. `initial` pre-computes initial estimates. `ic` comptes solution based on information criteria (BIC, AIC or AICc). `cv` computes solution based on cross-validation (cv). 
#' @param nlam number of \eqn{\lambda}'s to use in the regularization path.
#' @param lambdas user specified sequence of \eqn{\lambda} values for fitting. We recommend leaving this to NULL and letting function to self-select values.
#' @param min_frac the minimum value of the penalty parameter, as a fraction of the maximum value.
#' @param nfolds number of folds of the cv loop.
#' @param lambda_choice chose between `min` and `1se`. `min` computes solution that minimizes the cv error. `1se` computes solution such that the cv error is within 1 standard error of the minimum `min`.
#' @param ic_choice choose between `bic`, `aic` and `aicc`. `bic` computes solution that minimizes Bayesian information criterion. `aic` computes solution that minimizes Akaike information criterion. `aicc` omputes solution that minimizes Akaike information corrected criterion.
#' @param num_cores number of cores used to compute cv loop.
#' @param verbose flag to print information. 
#' @param thresh convergence threshold for change in beta. We recommend leaving this to NULL.
#' @param outer_thresh outer loop convergence threshold. We recommend leaving this to NULL.
#' @param inner_iter the maximum number of inner sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @param outer_iter the maximum number of outer sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @return Parameter estimates of panel data regression model under sg-LASSO penalty.
#' @author Jonas Striaukas
#' @examples 
#' \donttest{
#' # simulate DGP
#' set.seed(1)
#' t <- 21; n = 20; p = 100; size.groups = 4 
#' index <- ceiling(1:p / size.groups)
#' X <- matrix(rnorm(n * t * p), ncol = p, nrow = n*t)
#' beta <- c(5,4,3,2,1)
#' y <- X[,1:5] %*% beta + 5*rnorm(n*t)
#' entity_indices <- sort(rep(1:n,times=t-1))
#' panel_sgl(X = X, Z = NULL, y = y, index = index, 
#'   entity_indices = entity_indices, gamma_w = 1, 
#'   regress_choice = "fe", method_choice = "initial", 
#'   num_cores = 2, verbose = FALSE)
#' }
#' @export panel_sgl
panel_sgl <- function(X, Z=NULL, y, index, entity_indices, gamma_w=NULL, l1_factor=NULL, l21_factor=NULL, dummies_index=NULL, full_est=NULL, regress_choice=c("re","fe"), method_choice=c("ic","cv","initial"), nlam=100, lambdas=NULL,min_frac=NULL, nfolds=10, lambda_choice=c("min","1se"), ic_choice=c("bic","aic","aicc"),
                            num_cores = NULL, verbose=FALSE,thresh=NULL, outer_thresh=NULL, inner_iter=NULL, outer_iter=NULL){
  # check if input data has no na entries
  if(any(is.na(y))){stop("y has NA entries, check and rerun")}
  if(any(is.na(X))){stop("X has NA entries, check and rerun")}
  method <- match.arg(method_choice)
  reg <- match.arg(regress_choice)
  lambda_choice <- match.arg(lambda_choice)
  
  ic_choice <- match.arg(ic_choice)
  
  if (lambda_choice=="min") {which_lambda = as.double(0)}
  if (lambda_choice=="1se") {which_lambda = as.double(1)}
  
  
  ntp <- dim(X)
  nt <- ntp[1]
  p <- ntp[2]
  n <- max(entity_indices)
  t <- nt/n
  
  # check if t and n are the whole numbers:
  checkfun <- function(input,tol=.Machine$double.eps)
    min(abs(c(input%%1, input%%1-1))) < tol
  if(checkfun(t)==FALSE)
    stop("sample size is not a whole number. check the input matrix X")
  if(checkfun(n)==FALSE)
    stop("number of entities is not a whole number. check input vector entity_indices")
  
  # sort dummies
  if(is.null(Z)){
    if (reg=="re"){
      Z <- rep(1,times=nt)
    }
    if (reg=="fe"){ 
      Z <- kronecker(diag(n),rep(1,times=t))
    }
  }
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  y <- as.vector(y)
  # check Z
  if (reg=="re"){
    if (dim(Z)[2]>1)
      stop("wrong dummy matrix was inputed for Random Effects regression. Please check and re-run.")
  }
  if (reg=="fe"){
    if (dim(Z)[2]!=n)
      stop("wrong dummy matrix was inputed for Fixed Effects regression. Please check and re-run.")
  }
  nparams <- dim(Z)[2]
  # sort penalty terms
  if (reg=="re"){
    if (is.null(l1_factor)){
      l1_factor <- 0
    }
    if (is.null(l21_factor)){
      l21_factor <- 0
    }
    if (is.null(dummies_index)){
      dummies_index <- 1
    }
  }
  if (reg=="fe"){
    if (is.null(l1_factor)){
      l1_factor <- rep(0,times=n)
    }
    if (is.null(l21_factor)){
      l21_factor <- rep(0,times=n)
    }
    if (is.null(dummies_index)){
      dummies_index <- 1:n
    }
  }
  if (is.null(lambdas)){
    if(verbose)
      message(paste0("computing lambda sequence"))
    lambdas <- path_calc_panel(X, Z, y, index, gamma_w, l1_factor, l21_factor, dummies_index, nlam = nlam, min_frac = min_frac)
  }
  # compute entire solution
  if (is.null(full_est)){
    if(verbose)
      message(paste0("computing full path solution"))
    full_est <- sgl_fit(X, Z, y, index, lambdas, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter, outer_iter, thresh, outer_thresh)
  }
  if (method=="initial"){
    if(verbose)
      message(paste0("returning initial full path estimates and lambda sequence"))
    return(list(full_est=full_est,lambdas=lambdas))
  }
  if (method=="ic"){
    if(verbose) 
      message(paste0("computing solution using - ",ic_choice," - as information criterion"))
      alpha_path <- full_est$alphahat
      beta_path <- full_est$betahat
      if (reg=="fe"){
        if (dim(alpha_path)[2]!= dim(beta_path)[2])
          stop("alpha and beta have different number of path estimates. check lambda configuration.")
        if (nlam!=dim(alpha_path)[2] || nlam!=dim(beta_path)[2]){
          nlam <- dim(alpha_path)[2]
          if(verbose)
            message("nlam was reset to match the number of coefficient estimates in path solution. perhaps non-default nlam was used to fit initial model without re-setting the value for ic computation")
        }
      }
      if (reg=="re"){
        if (nlam!=dim(beta_path)[2]){
          nlam <- dim(beta_path)[2]
          if(verbose)
            message("nlam was reset to match the number of coefficient estimates in path solution. perhaps non-default nlam was used to fit initial model without re-setting the value for ic computation")
        }
      }
      crit <- numeric(nlam)
    if (reg=="re"){
      alpha_path <- matrix(alpha_path,nrow=1,ncol=nlam)
    } 
    for (i in 1:nlam){
      tmp <- NULL
      tmp$alpha <- alpha_path[,i]
      tmp$beta <- beta_path[,i]
      # compute df as |beta|_0
      beta_0_card <- sum(tmp$beta!=0)
      df <- beta_0_card+nparams
      # compute penalty
      if (ic_choice=="bic") {pen <- log(nt)/nt*df}
      if (ic_choice=="aic") {pen <- 2/nt*df}
      if (ic_choice=="aicc") {pen <- 2*df/(nt - df - 1)}
      yhat <- predict.panel_sgl(object = tmp, newX = X, newZ = Z, regress_choice = reg)$pred
      sigsqhat <- sum((y-mean(y))^2)/nt
      mse <- sum((y-yhat)^2)/nt
      crit[i] <- mse/sigsqhat + pen 
    }
    lambda_ic <- lambdas[which(crit==min(crit))]
    #take single solution in case many are optimum:
    if (length(lambda_ic)>1){lambda_ic <- min(lambda_ic)} 
    min_crit <- list(crit)
    beta <- full_est$betahat[,which(lambdas==lambda_ic)]
    if (reg=="fe"){
      alpha <-  full_est$alphahat[,which(lambdas==lambda_ic)]
    } 
    if (reg=="re"){
      alpha <-  full_est$alphahat[which(lambdas==lambda_ic)]
    } 
    return(list(alpha_path = alpha_path, beta_path = beta_path, lambdas = lambdas, alpha = alpha, beta = beta, lambda_ic = lambda_ic, min_crit = min_crit, full_est = full_est))
  }
  if (method=="cv") {
    if(verbose) {message(paste0("computing ", nfolds, "-fold cross-validation"))}
    # compute entity-specific folds
    foldsid <- numeric(nt)
    for (i in 1:n){
      idx <- which(entity_indices==i)
      i_nrow <- length(idx)
      foldsid[idx] <- cvfolds(as.double(nfolds), as.double(i_nrow))+1
    }
    #======== mc stuff ========# 
    if (is.null(num_cores)){
      if(.Platform$OS.type=="unix")
        num_cores <- parallel::detectCores()
      if(.Platform$OS.type=="windows")
        num_cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    }
    cl <- parallel::makeCluster(num_cores)
    registerDoSNOW(cl)
    if(verbose) {
      pb <- utils::txtProgressBar(max=nfolds, style=3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
    } else {
      opts <- NULL
    }
    #======== Main loop ========# 
    output <- foreach::foreach(k = 1:nfolds, .packages = c("midasml"), .options.snow=opts) %dopar% {
      # computing CV solutions
      which_val <- which(foldsid==k)
      which_train <- which(foldsid!=k)
      
      y_train <- y[which_train]
      X_train <- X[which_train, ]
      Z_train <- Z[which_train, ]
      
      est_new <- sgl_fit(X_train, Z_train, y_train, index, lambdas, gamma_w, l1_factor, l21_factor, dummies_index)
      est_new
    }
    parallel::stopCluster(cl)
    # initialize cv mean and sd storage
    cvm <- numeric(nlam)
    cvsd_fold <- matrix(0, nrow=nlam, ncol=nfolds)
    for (k in 1:nfolds){
      tmp <- output[[k]]
      est_new <- rbind(tmp$alphahat,tmp$betahat)
      which_val <- which(foldsid==k)
      which_train <- which(foldsid!=k)
      
      y_val <- y[which_val]
      X_val <- X[which_val, ]
      Z_val <- Z[which_val, ]
      ZX_val <- as.matrix(cbind(Z_val,X_val))
      tmp <- updatecvmsd(as.vector(cvm),  as.matrix(cvsd_fold), as.double(nlam), as.matrix(est_new), as.double(k), as.vector(y_val), as.matrix(ZX_val))
      cvm <- tmp$cvm
      cvsd <- tmp$cvsd_fold
    }
    cvm <- cvm/nfolds
    cvsd <- apply(cvsd,1,stats::sd) * sqrt(nfolds)
    lambda_cv <- getmin_cpp(lambdas, cvm, cvsd, which_lambda)
    min_crit <- list(cvm = cvm, cvsd = cvsd)
    alpha_path <- full_est$alphahat
    beta_path <- full_est$betahat
    if (reg=="fe")
      alpha <- full_est$alphahat[,which(lambdas==lambda_cv)]
    if (reg=="re")
      alpha <- full_est$alphahat[which(lambdas==lambda_cv)]
    
    beta <- full_est$betahat[,which(lambdas==lambda_cv)]
    
    return(list(alpha_path = alpha_path, beta_path = beta_path, lambdas = lambdas, alpha = alpha, beta = beta, lambda_cv = lambda_cv, min_crit = min_crit, full_est = full_est))
  } 
}

#' Linear sg-LASSO regression
#' 
#' @description 
#' Fits sg-LASSO least squares regression model. Options include cross-validation and information criteria for \eqn{\lambda} penalty parameter selection.
#' @details
#' \ifelse{html}{\out{The sequence of linear regression models implied by <code>lambdas</code> vector is fit by block coordinate-descent. The objective function is <br><br> <center> RSS(&alpha;,&beta;)/T + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where RSS(&alpha;,&beta;) is the least squares fit. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The sequence of linear regression models implied by \code{lambdas} vector is fit by block coordinate-descent. The objective function is \deqn{RSS(\alpha,\beta)/T + 2\lambda * \Omega_\gamma(\beta),} where \eqn{RSS(\alpha,\beta)} is the least squares fit. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage 
#' reg_sgl(X, y, index, gamma_w = NULL, full_est = NULL, 
#'   method_choice = c("ic", "cv", "initial"), nlam = 100, lambdas = NULL, 
#'   min_frac = NULL, nfolds = 10, lambda_choice = c("min", "1se"), 
#'   ic_choice = c("bic", "aic", "aicc"), num_cores = NULL, verbose = FALSE, 
#'   thresh = NULL, outer_thresh = NULL, inner_iter = NULL, outer_iter = NULL)
#' @param X T by p data matrix, where t and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 vector of outcome.
#' @param index p by 1 vector indicating group membership of each covariate.
#' @param gamma_w sg-LASSO mixing parameter. \code{gamma_w = 1} is LASSO and \code{gamma_w = 0} group LASSO.
#' @param full_est pre-estimated parameters based on full sample and \code{regress_choice} for a sequence of \eqn{\lambda}'s.
#' @param method_choice choose between `initial`, `ic` and `cv`. `initial` pre-computes initial estimates. `ic` comptes solution based on information criteria (BIC, AIC or AICc). `cv` computes solution based on cross-validation (cv). 
#' @param nlam number of \eqn{\lambda}'s to use in the regularization path.
#' @param lambdas user specified sequence of \eqn{\lambda} values for fitting. We recommend leaving this to NULL and letting function to self-select values.
#' @param min_frac the minimum value of the penalty parameter, as a fraction of the maximum value.
#' @param nfolds number of folds of the cv loop.
#' @param lambda_choice chose between `min` and `1se`. `min` computes solution that minimizes the cv error. `1se` computes solution such that the cv error is within 1 standard error of the minimum `min`.
#' @param ic_choice choose between `bic`, `aic` and `aicc`. `bic` computes solution that minimizes Bayesian information criterion. `aic` computes solution that minimizes Akaike information criterion. `aicc` omputes solution that minimizes Akaike information corrected criterion.
#' @param num_cores number of cores used to compute cv loop.
#' @param verbose flag to print information. 
#' @param thresh convergence threshold for change in beta. We recommend leaving this to NULL.
#' @param outer_thresh outer loop convergence threshold. We recommend leaving this to NULL.
#' @param inner_iter the maximum number of inner sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @param outer_iter the maximum number of outer sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @return Parameter estimates of linear regression model under sg-LASSO penalty.
#' @author Jonas Striaukas
#' @examples 
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' index = 1:20
#' reg_sgl(X = x, y = y, index = index, gamma_w = 1, method_choice = "initial", 
#'   num_cores = 2, verbose = FALSE)
#' @export reg_sgl
reg_sgl <- function(X, y, index, gamma_w=NULL, full_est=NULL, method_choice=c("ic","cv","initial"), nlam=100, lambdas=NULL,min_frac=NULL, nfolds=10, lambda_choice=c("min","1se"), ic_choice=c("bic","aic","aicc"),
                        num_cores = NULL, verbose=FALSE,thresh=NULL, outer_thresh=NULL, inner_iter=NULL, outer_iter=NULL){
  # check if input data has no na entries
  if(any(is.na(y)))
    stop("y has NA entries, check and rerun")
  if(any(is.na(X)))
    stop("X has NA entries, check and rerun")
  # get settings
  method <- match.arg(method_choice)
  lambda_choice <- match.arg(lambda_choice)
  ic_choice <- match.arg(ic_choice)
  
  if (lambda_choice=="min")
    which_lambda = as.double(0)
  if (lambda_choice=="1se")
    which_lambda = as.double(1)
  
  tp <- dim(X)
  t <- tp[1]
  p <- tp[2]
  
  # intercept dummy
  Z <- as.matrix(rep(1,times=length(y)))
  # zero penalty weight for the intercept
  l1_factor <- 0
  l21_factor <- 0 
  dummies_index <- 1 # group index for the intercept
  X <- as.matrix(X)
  y <- as.vector(y)
  if (is.null(lambdas)){
    if(verbose)
      message(paste0("computing lambda sequence"))
    lambdas <- path_calc(X, Z, y, index, gamma_w, l1_factor, l21_factor, dummies_index)
  }
  # compute entire solution
  if (is.null(full_est)){
    if(verbose)
      message(paste0("computing full path solution"))
    full_est <- sgl_fit(X, Z, y, index, lambdas, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter, outer_iter, thresh, outer_thresh)
  }
  if (method=="initial"){
    if(verbose) {message(paste0("returning initial full path estimates and lambda sequence"))}
    return(list(full_est = full_est, lambdas = lambdas))
  }
  if (method=="ic"){
    if(verbose) {message(paste0("computing solution using - ",ic_choice," - as information criterion"))}
    crit <- numeric(nlam)
    alpha_path <- full_est$alphahat
    beta_path <- full_est$betahat
    for (i in 1:nlam){
      tmp <- NULL
      tmp$alpha <- alpha_path[i]
      tmp$beta <- beta_path[,i]
      # compute df as |beta|_0
      beta_0_card <- sum(tmp$beta != 0)
      df <- beta_0_card
      # compute penalty
      if (ic_choice=="bic") {pen <- log(t)/t*df}
      if (ic_choice=="aic") {pen <- 2/t*df}
      if (ic_choice=="aicc") {pen <- 2*df/(t - df - 1)}
      yhat <- predict.reg_sgl(object = tmp, newX = X)$pred
      sigsqhat <- sum((y-mean(y))^2)/t
      mse <- sum((y-yhat)^2)/t
      crit[i] <- mse/sigsqhat + pen 
    }
    lambda_ic <- lambdas[which(crit==min(crit))]
    #take single solution in case many are optimum:
    if (length(lambda_ic)>1)
      lambda_ic <- min(lambda_ic)
    min_crit <- list(crit)
    alpha <- full_est$alphahat[which(lambdas==lambda_ic)]
    beta <- full_est$betahat[,which(lambdas==lambda_ic)]
    return(list(alpha_path = alpha_path, beta_path = beta_path, lambdas = lambdas, alpha = alpha, beta = beta, lambda_ic = lambda_ic, min_crit = min_crit, full_est = full_est))
  }
  if (method=="cv") {
    if(verbose)
      message(paste0("computing ", nfolds, "-fold cross-validation"))
    # compute folds
    foldsid <- cvfolds(as.double(nfolds), as.double(t))+1

    #======== mc stuff ========# 
    if (is.null(num_cores)){
      if(.Platform$OS.type=="unix")
        num_cores <- parallel::detectCores()
      if(.Platform$OS.type=="windows")
        num_cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    }
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    if(verbose) {
      pb <- utils::txtProgressBar(max=nfolds, style=3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
    } else {
      opts <- NULL
    }
    #======== Main loop ========# 
    output <- foreach::foreach(k = 1:nfolds, .packages = c("midasml"), .options.snow=opts) %dopar% {
        # computing CV solutions
        which_val <- which(foldsid==k)
        which_train <- which(foldsid!=k)

        y_train <- y[which_train]
        X_train <- X[which_train, ]
        Z_train <- Z[which_train, ]
        
        est_new <- sgl_fit(X_train, Z_train, y_train, index, lambdas, gamma_w, l1_factor, l21_factor, dummies_index)
        est_new
    }
    parallel::stopCluster(cl)
    # initialize cv mean and sd storage
    cvm <- numeric(nlam)
    cvsd_fold <- matrix(0, nrow=nlam, ncol=nfolds)
    for (k in 1:nfolds){
      tmp <- output[[k]]
      est_new <- rbind(tmp$alphahat,tmp$betahat)
      which_val <- which(foldsid==k)
      which_train <- which(foldsid!=k)
      
      y_val <- y[which_val]
      X_val <- X[which_val, ]
      Z_val <- Z[which_val, ]
      ZX_val <- as.matrix(cbind(Z_val,X_val))
      tmp <- updatecvmsd(as.vector(cvm),  as.matrix(cvsd_fold), as.double(nlam), as.matrix(est_new), as.double(k), as.vector(y_val), as.matrix(ZX_val))
      cvm <- tmp$cvm
      cvsd <- tmp$cvsd_fold
    }
    cvm <- cvm/nfolds
    cvsd <- apply(cvsd,1,sd) * sqrt(nfolds)
    lambda_cv <- getmin_cpp(lambdas, cvm, cvsd, which_lambda)
    min_crit <- list(cvm = cvm, cvsd = cvsd)
    alpha_path <- full_est$alphahat
    beta_path <- full_est$betahat
    alpha <- full_est$alphahat[which(lambdas==lambda_cv)]
    beta <- full_est$betahat[,which(lambdas==lambda_cv)]
    return(list(alpha_path = alpha_path, beta_path = beta_path, lambdas = lambdas, alpha = alpha, beta = beta, lambda_cv = lambda_cv, min_crit = min_crit, full_est = full_est))
  } 
}

#' sg-LASSO regression
#' @description 
#' Fits mse sg-LASSO regression model. 
#' @details
#' \ifelse{html}{\out{The sequence of linear regression models implied by <code>lambdas</code> vector is fit by block coordinate-descent. The objective function is <br><br> <center> RSS(&alpha;,&beta;)/T + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where RSS(&alpha;,&beta;) is the least squares fit. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The sequence of linear regression models implied by \code{lambdas} vector is fit by block coordinate-descent. The objective function is \deqn{RSS(\alpha,\beta)/T + 2\lambda * \Omega_\gamma(\beta),} where \eqn{RSS(\alpha,\beta)} is the least squares fit. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage
#' sgl_fit(X, Z, y, index, lambdas, gamma_w = NULL, l1_factor = NULL, l21_factor = NULL,
#'   dummies_index = NULL, inner_iter = NULL, outer_iter = NULL, thresh = NULL, 
#'   outer_thresh = NULL)
#' @param X T by p data matrix, where t and p respectively denote the sample size and the number of regressors.
#' @param Z dummies matrix.
#' @param y T by 1 vector of outcome.
#' @param index p by 1 vector indicating group membership of each covariate.
#' @param lambdas user specified sequence of \eqn{\lambda} values for fitting. We recommend leaving this to NULL and letting function to self-select values.
#' @param gamma_w sg-LASSO mixing parameter. \code{gamma_w = 1} is LASSO and \code{gamma_w = 0} group LASSO.
#' @param l1_factor \ifelse{html}{\out{&#8467;<sub>1</sub> norm penalty factor for random or fixed effects (default value is zero which means &alpha; is left unpenalized in &#8467;<sub>1</sub> norm).}}{\eqn{\ell_1} norm penalty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1} norm).}     
#' @param l21_factor \ifelse{html}{\out{&#8467;<sub>2,1</sub> norm penalty factor for random or fixed effects (default value is zero which means &alpha; is left unpenalized in &#8467;<sub>2,1</sub> norm).}}{\eqn{\ell_1} norm penalty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1} norm).}     
#' @param dummies_index vector indicating group membership of \eqn{\alpha} (default - no grouping).
#' @param thresh convergence threshold for change in beta. We recommend leaving this to NULL.
#' @param outer_thresh outer loop convergence threshold. We recommend leaving this to NULL.
#' @param inner_iter the maximum number of inner sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @param outer_iter the maximum number of outer sg-LASSO loop iterations. We recommend leaving this to NULL.
#' @return sg-LASSO regression fitted coefficients.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' index = 1:20
#' Z <- as.matrix(rep(1,times=length(y)))
#' sgl_fit(X = x, Z = Z, y = y, index = index, lambdas = c(1,2,3), gamma_w = 1)
#' @export sgl_fit
sgl_fit <- function(X, Z, y, index, lambdas, gamma_w=NULL, l1_factor=NULL, l21_factor=NULL, dummies_index=NULL, inner_iter=NULL, outer_iter=NULL, thresh=NULL, outer_thresh=NULL){
  nlam <- length(lambdas)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  y <- as.vector(y)
  num_groups = max(index)
  if(is.null(lambdas)) 
    stop("lambda sequence was not specified")
  if(is.null(l1_factor)) 
    l1_factor <- rep(0,dim(Z)[2])
  if(is.null(l21_factor)) 
    l21_factor <- rep(0,dim(Z)[2])
  if(is.null(dummies_index)) 
    dummies_index <- 1
  
  if(is.null(gamma_w)){
    message("sg-lasso relative weight was not set. setting to default value 0.8")
    gamma_w <- 0.8
  }
  
  if(is.null(inner_iter))
    inner_iter <- 1e2
  
  if(is.null(outer_iter))
    outer_iter <- 1e4
  
  if(is.null(thresh))
    thresh <- 1e-2
  
  if(is.null(outer_thresh))
    outer_thresh <- 1e-4
  
  gamma_solver <- 0.8
  step <- 1
  reset <- 10
  
  zparam <- dim(Z)[2]
  xparam <- dim(X)[2]
  dummies <- 1
  
  #------------------ sgl fit -----------#
  fit <- cpp_sgl_fitpath(as.matrix(X),as.matrix(Z),as.vector(y), as.vector(index), as.double(dummies), 
                  as.vector(l1_factor),as.vector(l21_factor),as.vector(dummies_index),
                  as.vector(lambdas), as.double(gamma_w), as.integer(inner_iter), as.integer(outer_iter), as.double(thresh), as.double(outer_thresh), 
                  as.double(gamma_solver), as.double(step), reset = as.integer(reset))
  alphahat <- fit[1:zparam,]
  betahat <- fit[-c(1:zparam),]
  return(list(alphahat=alphahat,betahat=betahat))
}

#' Computes prediction for the sg-LASSO panel regression model
#' 
#' @param object fit object from \code{panel_sgl}.
#' @param newX matrix of out-of-sample covariate observations.
#' @param newZ optional matrix of dummies for panel data model.
#' @param regress_choice choose between `re` and `fe`. Must be consistent with \code{object}.
#' @param ... currently ignored optional parameters. 
#' @return a list of these variables:
#' @return pred - overall prediction.
#' @return predZ - dummies prediction.
#' @return predX - covariates prediction.
#' @examples 
#' set.seed(1)
#' t <- 21; n = 20; p = 100; size.groups = 4 
#' index <- ceiling(1:p / size.groups)
#' X <- matrix(rnorm(n * t * p), ncol = p, nrow = n*t)
#' beta <- c(5,4,3,2,1)
#' y <- X[,1:5] %*% beta + 5*rnorm(n*t)
#' Z <- kronecker(diag(n), rep(1, times = t))
#' entity_indices <- sort(rep(1:n,times=t-1))
#' fit <- panel_sgl(X = X, Z = Z, y = y, index = index, 
#'          entity_indices = entity_indices, gamma_w = 1, 
#'          regress_choice = "fe", method_choice = "ic", 
#'          num_cores = 2, verbose = FALSE)
#' predict.panel_sgl(object = fit, newX = X, newZ = Z, regress_choice = "fe")$pred
#' @author Jonas Striaukas
#' @method predict panel_sgl
#' @rdname predict.panel_sgl
#' @export predict.panel_sgl
predict.panel_sgl <- function(object, newX, newZ=NULL, regress_choice=c("re","fe"), ...){
  newX <- as.matrix(newX)
  n <- dim(newX)[1]
  reg <- match.arg(regress_choice)
  
  alpha <- object$alpha
  beta <- object$beta
  if (is.null(newZ)){
    if (reg=="re"){
      Z <- matrix(rep(1,times=n),nrow=n,ncol=1)
    }
    if (reg=="fe"){ 
      Z <- diag(n)
    }
  } else {
    Z <- newZ
  }
  predZ <- as.numeric(Z%*%alpha)
  predX <- as.numeric(newX%*%beta)
  pred <- as.numeric(predZ+predX)
  return(list(pred=pred,predZ=predZ,predX=predX))
}

#' Computes prediction for the sg-LASSO linear regression
#' 
#' @param object fit object from sglassofit.
#' @param newX matrix of out-of-sample covariate observations.
#' @param ... currently ignored optional parameters. 
#' @return a list of these variables:
#' @return pred - overall prediction.
#' @return predZ - intercept prediction.
#' @return predX - covariates prediction.
#' @author Jonas Striaukas
#' @examples 
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' index <- 1:20
#' fit <- reg_sgl(X = x, y = y, index = index, gamma_w = 1, method_choice = "ic", 
#'        num_cores = 2, verbose = FALSE)
#' predict.reg_sgl(object = fit, newX = x)
#' @method predict reg_sgl
#' @rdname predict.reg_sgl
#' @export predict.reg_sgl
predict.reg_sgl <- function(object, newX, ...){
  alpha <- object$alpha
  beta <- object$beta
  predZ <- alpha
  predX <- as.numeric(newX%*%beta)
  pred <- as.numeric(predZ+predX)
  return(list(pred=pred,predZ=predZ,predX=predX))
}

#' Computes a sequence of lambda parameters for the mse sg-LASSO regression
#' 
#' @param X T by p data matrix, where t and p respectively denote the sample size and the number of regressors.
#' @param Z dummies matrix left for the interecept. 
#' @param y T by 1 vector of outcome.
#' @param index p by 1 vector indicating group membership of each covariate.
#' @param gamma_w sg-LASSO mixing parameter. \code{gamma_w = 1} is LASSO and \code{gamma_w = 0} group LASSO.
#' @param l1_factor penalty term on the intercept. Should be left unpenalized, i.e set to zero (default).
#' @param l21_factor penalty term on the intercept. Should be left unpenalized, i.e set to zero (default).
#' @param dummies_index vector indicating dummies group membership. Should be left unspecified, i.e set to one (default).
#' @param nlam number of \eqn{\lambda}'s to use in the regularization path.
#' @param min_frac the minimum value of the penalty parameter, as a fraction of the maximum value.
#' @return lambdas sequence of \eqn{\lambda} values for fitting.
#' @export path_calc
#' @keywords internal
path_calc <- function(X, Z, y, index, gamma_w=NULL, l1_factor=0, l21_factor=0, dummies_index=1, nlam=NULL, min_frac=NULL){
  if(is.null(gamma_w)){
    message("sg-lasso relative weight was not set. setting to default value 0.8")
    gamma_w <- 0.8
  }
  if(is.null(nlam))
    nlam <- 100
  
  if (is.null(min_frac)){
    k <- dim(X)[2]
    t <- dim(X)[1]
    if (k>t){
      min_frac <- 0.0001
    } else {
      min_frac <- 0.001
    }
  }
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  y <- as.vector(y)
  n <- dim(X)[1]
  max_lam_tmp <- max(t(cbind(Z,X))%*%y)/n
  
  # we just need approximate solutions so coordinate-descent parameters are set to be loose
  inner_iter <- 1e3
  outer_iter <- 1e3
  thresh <- 1e-3
  outer_thresh <- 1e-3
  
  # fit the initial solution
  path_est <- sgl_fit(X, Z, y, index, max_lam_tmp, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
  tmp_lambda <- max_lam_tmp
  notfound <- TRUE
  k <- 0
  if(any(path_est!=0)){
    factor <- 1.01
    while (notfound){
      k <- k + 1
      tmp_lambda <- tmp_lambda*factor
      tmp_fit <- sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
      notfound <- any(tmp_fit!=0)
      if (k > 1e3){
        notfound <- FALSE
        tmp_lambda <- max_lam_tmp
        warning("more than 1e3 iterations were used to fine-tune lambdas sequence. initial value set to max(Xty)/n")
      }
    }
  } else {
    factor <- 0.99
    while (notfound){
      k <- k + 1
      tmp_lambda <- tmp_lambda*factor
      tmp_fit <- sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
      notfound <- all(tmp_fit==0)
      if (k > 1e3){
        notfound <- FALSE
        tmp_lambda <- max_lam_tmp
        warning("more than 1e3 iterations were used to fine-tune lambdas sequence. initial value set to max(Xty)/n")
      }
    }
    tmp_lambda <- tmp_lambda/0.99
  }
  max_lam <- tmp_lambda
  min_lam <- min_frac*max_lam
  log_min_lam <- log(min_lam)
  log_max_lam <- log(max_lam)
  seq_loglambdas <- seq(log_max_lam,log_min_lam,length.out = nlam)
  lambdas <- exp(seq_loglambdas)
  lambdas
}

#' Computes a sequence of lambda parameters for the sg-LASSO panel regression
#' 
#' @param X nT by p data matrix, where n, t and p respectively denote the number of individuals, sample size and the number of regressors.
#' @param Z dummies matrix for random effects or fixed effects panel data model. If left unspecified, it is computed based on \code{regress_choice} choice.
#' @param y nT by 1 vector of outcome.
#' @param index p by 1 vector indicating group membership of each covariate.
#' @param gamma_w sg-LASSO mixing parameter. \code{gamma_w = 1} is LASSO and \code{gamma_w = 0} group LASSO.
#' @param l1_factor \eqn{\ell_1} norm peanlty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1} norm).
#' @param l21_factor \eqn{\ell_{1,2}} norm peanlty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_{1,2}} norm).
#' @param dummies_index vector indicating group membership of \eqn{\alpha}.
#' @param nlam number of \eqn{\lambda}'s to use in the regularization path.
#' @param min_frac the minimum value of the penalty parameter, as a fraction of the maximum value.
#' @return lambdas sequence of \eqn{\lambda} values for fitting.
#' @export path_calc_panel
#' @keywords internal
path_calc_panel <- function(X, Z, y, index, gamma_w=NULL, l1_factor=NULL, l21_factor=NULL, dummies_index=NULL, nlam=NULL, min_frac=NULL){
  if(is.null(gamma_w)){
    message("sg-lasso relative weight was not set. setting to default value 0.8")
    gamma_w <- 0.8
  }
  if(is.null(nlam))
    nlam <- 100
  
  if (is.null(min_frac)){
    k <- dim(X)[2]
    nt <- dim(X)[1]
    t <- nt/dim(Z)[2]
    if (k>t){
      min_frac <- 0.0001
    } else {
      min_frac <- 0.001
    }
  }
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  y <- as.vector(y)
  n <- dim(X)[1]
  max_lam_tmp <- max(t(cbind(Z,X))%*%y)/n
  
  # we just need approximate solutions so coordinate-descent parameters are set to be loose
  inner_iter <- 1e1
  outer_iter <- 1e1
  thresh <- 1e-1
  outer_thresh <- 1e-1
  
  # fit the initial solution
  path_est <- sgl_fit(X, Z, y, index, max_lam_tmp, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
  tmp_lambda <- max_lam_tmp
  notfound <- TRUE
  if(any(path_est!=0)){
    factor <- 1.05
    while (notfound){
      tmp_lambda <- tmp_lambda*factor
      tmp_fit <- sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
      notfound <- any(tmp_fit!=0)
    }
  } else {
    factor <- 0.95
    while (notfound){
      tmp_lambda <- tmp_lambda*factor
      tmp_fit <- sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor, dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, thresh=thresh, outer_thresh=outer_thresh)$betahat
      notfound <- all(tmp_fit==0)
    }
    tmp_lambda <- tmp_lambda/0.95
  }
  max_lam <- tmp_lambda
  min_lam <- min_frac*max_lam
  log_min_lam <- log(min_lam)
  log_max_lam <- log(max_lam)
  seq_loglambdas <- seq(log_max_lam,log_min_lam,length.out = nlam)
  lambdas <- exp(seq_loglambdas)
  lambdas
}

#' Updates the cross-validation error estimates
#' 
#' @param cvm current cvm estimates.
#' @param cvsd_fold current cvsd_fold estimates. 
#' @param nlam number of \eqn{\lambda}'s to use in the regularization path.
#' @param est_new p by nlam matrix sequence of estimates, where p is the total dimension of the data matrix and nlam is the number of \eqn{\lambda}'s to use in the regularization path.
#' @param k fold index. 
#' @param y_val response used in validation.  
#' @param ZX_val y_val data matrix used in validation.  
#' @return updated cross-validation error estimates.
#' @export updatecvmsd
#' @keywords internal
updatecvmsd <- function(cvm,  cvsd_fold, nlam, est_new, k, y_val, ZX_val){
  t <- length(y_val)
  for (d in seq(nlam)){
    beta_d <- est_new[,d]
    cvsd_fold[d,k] = sum((y_val - ZX_val*beta_d)^2)/t
    cvm[d] = cvm[d] + sum((y_val - ZX_val*beta_d)^2)/t
  }
  return(list(cvm=cvm,cvsd_fold=cvsd_fold))
}

