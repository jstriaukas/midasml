# Model: ADL-MIDAS (FADL-MIDAS) regression 
# Estimation: MIDAS-NLS based on analytic derivatives (gradient and Hessian)
# Polynomial: Beta, 3 parametrizations, refer to "poly_spec"
# Paper, analytic derivatives:  
# - Kostrov (2020) Estimating MIDAS regressions via MIDAS-NLS with revised optimization. Working paper. 
# Optimizer: "nlminb", with constraints. 
optim_ardl_beta <- function(y, z, x, poly_spec = 0, nbtrials = 100){
  n <- length(y)
  k <- dim(z)[2]
  if (is.null(k)){
    z <- matrix(z, nrow=length(z))
    k <- 1
  }
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)

  # store data into a list:
  dataList <- list (y=y, z=z, x=x, poly_spec=poly_spec)
  # get parameter constraints:
  constr <- get_constr_beta(poly_spec, k)
  #-------------------- main optimization --------------------#
  invisible(capture.output(opt <- multiStartoptim(objectivefn = estimate_ardl_beta, 
                                        gradient = gradient_ardl_beta,
                                        data = dataList,
                                        hessian = hessian_ardl_beta, 
                                        lower = constr$lw_b, upper = constr$up_b, 
                                        method = "nlminb", 
                                        nbtrials = nbtrials,
                                        typerunif = "runifbase",
                                        control = list(eval.max=1000, iter.max=1000,
                                        rel.tol=1e-12, x.tol=1.5e-10))))
  #-------------------- xxxxxxxxxxxxxxxx --------------------#
  # back-out parameters:
  k1 <- opt$par[1]
  k2 <- opt$par[2]
  k3 <- opt$par[3]
  beta <- opt$par[4]
  c <- opt$par[5]
  rho <- opt$par[6:length(opt$par)]

  
  # sort the output:
  other_x <- NULL
  for (i in 1:k) 
    other_x <- c(other_x, paste0("other-",i))
  
  int <- "(Intercept)"
  if (poly_spec == 0) {
    coef <- matrix(c(c, rho, beta, k1, k2, k3), nrow = 1, ncol = k + 5)
    rownames(coef) <- ""
    colnames(coef) <- c(int,other_x,"beta","k1","k2", "k3")
  } else if (poly_spec == 1) {
    coef <- matrix(c(c, rho, beta, k2, k3), nrow = 1, ncol = k + 4)
    rownames(coef) <- ""
    colnames(coef) <- c(int,other_x,"beta","k2", "k3")
  } else if (poly_spec == 2) {
    coef <- matrix(c(c, rho, beta, k1, k2), nrow = 1, ncol = k + 4)
    rownames(coef) <- ""
    colnames(coef) <- c(int,other_x,"beta","k1", "k2")
  }
  class(coef) <- "midas.beta"
  return(coef)
}

get_constr_beta <- function(poly_spec, k){
  # define upper and lower constraints
  up_theta1 <- 70
  up_theta2 <- 70
  up_theta3 <- 70
  up_beta <- 40
  up_rho <- 40
  
  low_theta1 <- 0
  low_theta2 <- 0
  low_theta3 <- 0
  low_beta <- -40
  low_rho <- -40
  
  # based on beta specification, construct vectors of constraints
  if (poly_spec==0){
    up_b = c(up_theta1, up_theta2, up_theta3, up_beta, rep(up_rho, k+1))
    lw_b = c(low_theta1, low_theta2, low_theta3, low_beta, rep(low_rho, k+1))  
  } else if (poly_spec==1){
    up_b = c(1, up_theta2, up_theta3, up_beta, rep(up_rho, k+1))
    lw_b = c(1, low_theta2, low_theta3, low_beta, rep(low_rho, k+1))  
  } else if (poly_spec==2){
    up_b = c(up_theta1, up_theta2, 0, up_beta, rep(up_rho, k+1))
    lw_b = c(low_theta1, low_theta2, 0, low_beta, rep(low_rho, k+1))  
  } else if (poly_spec==3){
    up_b = c(1, up_theta2, 0 , up_beta, rep(up_rho, k+1))
    lw_b = c(1, low_theta2, 0, low_beta, rep(low_rho, k+1)) 
  }
  return(list(up_b = up_b, lw_b = lw_b)) 
}

estimate_ardl_beta <- function(args, dataList) {
  
  y = dataList$y
  z = dataList$z
  x = dataList$x
  
  p = dim(x)[2]
  ii = matrix(1, p,1)     
  xx = matrix(c(1:p)/(p+1))  
  
  poly_spec = dataList$poly_spec
  
  
  if (poly_spec==0){
    theta1 = args[1]
    theta3 = args[3]
  } else if (poly_spec==1){
    theta1 = 1
    theta3 = args[3]     
  } else if (poly_spec==2){
    theta1 = args[1]
    theta3 = 0 
  } else if (poly_spec==3){
    theta1 = 1
    theta3 = 0    
  }
  
  theta2 = args[2]
  beta = args[4]
  rho = args[5:length(args)] 
  
  m = gamma(theta1+theta2)/(gamma(theta1)*gamma(theta2))
  weights = xx^(theta1-1)*(ii-xx)^(theta2-1)*m + theta3 
  #browser()
  mse = t(y  - z %*% rho -  beta * x %*% weights)%*%(y  - z %*% rho -  beta * x %*% weights)  
  
  return(as.numeric(mse))
  
}

gradient_ardl_beta <- function(args, dataList) {
  y = dataList$y
  z = dataList$z
  x = dataList$x
  
  p = dim(x)[2]
  ii = matrix(1, p,1)     
  xx = matrix(c(1:p)/(p+1))  
  
  poly_spec = dataList$poly_spec
  
  if (poly_spec==0){
    theta1 = args[1]
    theta3 = args[3]
  } else if (poly_spec==1){
    theta1 = 1
    theta3 = args[3]     
  } else if (poly_spec==2){
    theta1 = args[1]
    theta3 = 0 
  } else if (poly_spec==3){
    theta1 = 1
    theta3 = 0    
  }
  
  theta2 = args[2]
  beta = args[4]
  rho = args[5:length(args)] 
  m = gamma(theta1+theta2)/(gamma(theta1)*gamma(theta2))
  weights = xx^(theta1-1)*(ii-xx)^(theta2-1)*m + theta3 
  
  nabla_G = 2*(y- beta * x %*%  weights - z%*%rho )
  
  T_star_w_C =  -beta * t(x)
  T_star_beta_C =  -t(weights) %*% t(x)
  T_star_rho_C =  - t(z)
  
  weights_tilde = xx^(theta1-1)*(ii-xx)^(theta2-1)*m
  
  log_derivative_theta2 = log(ii-xx) + digamma(theta1+theta2) - digamma(theta2)
  T_star_theta2_W  = t(log_derivative_theta2*weights_tilde)
  part2 = T_star_theta2_W %*% T_star_w_C %*% nabla_G  # theta2
  
  if (poly_spec==0){
    log_derivative_theta1 = log(xx) + digamma(theta1+theta2) - digamma(theta1)
    T_star_theta1_W  = t(log_derivative_theta1*weights_tilde)
    T_star_theta3_W  = t(ii)
    part1 = T_star_theta1_W %*% T_star_w_C %*% nabla_G   # theta1
    part3 = T_star_theta3_W %*% T_star_w_C %*% nabla_G   # theta3
  } else if (poly_spec==1){
    T_star_theta3_W  = t(ii)
    part1 = 0   # theta1
    part3 = T_star_theta3_W %*% T_star_w_C %*% nabla_G   # theta3
  } else if (poly_spec==2){
    log_derivative_theta1 = log(xx) + digamma(theta1+theta2) - digamma(theta1)
    T_star_theta1_W  = t(log_derivative_theta1*weights_tilde)
    part1 = T_star_theta1_W %*% T_star_w_C %*% nabla_G   # theta1
    part3 = 0
  } else if (poly_spec==3){
    part1 = 0;   # theta1
    part3 = 0;   # theta3
  }
  
  part4 = T_star_beta_C  %*% nabla_G   # beta
  part5 = T_star_rho_C  %*% nabla_G    # rho
  grad =  c(part1, part2, part3, part4, part5)
  
  return(grad)  
}

hessian_ardl_beta <- function(args, dataList) {
  
  y = dataList$y
  z = dataList$z
  x = dataList$x
  
  p = dim(x)[2]
  ii = matrix(1, p,1)     
  xx = matrix(c(1:p)/(p+1)) 
  
  poly_spec = dataList$poly_spec
  
  
  if (poly_spec==0){
    theta1 = args[1]
    theta3 = args[3]
  } else if (poly_spec==1){
    theta1 = 1
    theta3 = args[3]     
  } else if (poly_spec==2){
    theta1 = args[1]
    theta3 = 0 
  } else if (poly_spec==3){
    theta1 = 1
    theta3 = 0    
  }
  
  theta2 = args[2]
  beta = args[4]
  rho = args[5:length(args)] 
  k = matrix(0, 1, length(rho))
  
  m = gamma(theta1+theta2)/(gamma(theta1)*gamma(theta2))
  weights = xx^(theta1-1)*(ii-xx)^(theta2-1)*m + theta3 
  
  ## line 1
  E1 = 2*(y - beta * x %*%  weights - z%*%rho )
  F1 =  -beta * t(x)
  D1 = log(xx) + digamma(theta1+theta2) - digamma(theta1)
  C1 = xx^(theta1-1)*(ii-xx)^(theta2-1)*m
  ### element 1
  dC1_theta1 = C1 * D1
  dD1_theta1 = trigamma(theta1+theta2 ) - trigamma(theta1)
  dB1_theta1 = F1 %*% (-2*beta*x%*%dC1_theta1)
  ### element 2
  D2 = log(ii-xx) + digamma(theta1+theta2) - digamma(theta2)
  dC1_theta2 = C1*D2
  dD1_theta2 = trigamma(theta1+theta2)
  dB1_theta2 = F1%*%(-2*beta*x%*%dC1_theta2)
  ### element 3
  dE1_theta3 = -2* beta*x %*% ii
  ### element 4
  dF1_beta = -t(x)
  dE1_beta = -2*x%*%weights
  ### element 5
  dE1_rho = -2*z
  
  ## line 2
  ### element 2
  dC2_theta2 = C1 * D2
  dD2_theta2 =  trigamma(theta1+theta2 ) - trigamma(theta2 )
  ### element 3
  ### element 4
  ### element 5
  H22 =  t(dC2_theta2 * D2 + C1 * dD2_theta2)%*%(F1%*%E1)  +  t(C1 *D2 )%*%dB1_theta2
  H23 = t(C1 * D2)%*%F1%*%dE1_theta3 
  H24 = t(C1*D2) %*% (dF1_beta%*%E1+F1%*%dE1_beta)
  H25 = t(C1 * D2)%*%F1%*%dE1_rho
  
  ## line 3
  A3 = t(ii) %*% F1
  ### element 3
  ### element 4
  dA3_beta = t(ii) %*% dF1_beta
  ### element 5
  
  if (poly_spec==0){
    H11 = t(dC1_theta1 * D1 + C1 * dD1_theta1)%*%(F1%*%E1)  +  t(C1 *D1 )%*%dB1_theta1
    H12 =  t(dC1_theta2 * D1 + C1 * dD1_theta2)%*%(F1%*%E1)  +  t(C1 *D1 )%*%dB1_theta2
    H13 = t(C1 * D1)%*%F1%*%dE1_theta3 
    H14 = t(C1 *D1) %*% (dF1_beta%*%E1+F1%*%dE1_beta)
    H15 = t(C1 * D1)%*%F1%*%dE1_rho
    H33 = A3 %*% dE1_theta3
    H34 = dA3_beta %*% E1 + A3 %*% dE1_beta
    H35 = A3 %*% dE1_rho 
  } else if (poly_spec==1){
    H11 = 0
    H12 = 0
    H13 = 0 
    H14 = 0
    H15 = k
    H33 = A3 %*% dE1_theta3
    H34 = dA3_beta %*% E1 + A3 %*% dE1_beta
    H35 = A3 %*% dE1_rho    
  } else if (poly_spec==2){
    H11 = t(dC1_theta1 * D1 + C1 * dD1_theta1)%*%(F1%*%E1)  +  t(C1 *D1 )%*%dB1_theta1
    H12 =  t(dC1_theta2 * D1 + C1 * dD1_theta2)%*%(F1%*%E1)  +  t(C1 *D1 )%*%dB1_theta2
    H13 = 0  
    H14 = t(C1 *D1) %*% (dF1_beta%*%E1+F1%*%dE1_beta)
    H15 = t(C1 * D1)%*%F1%*%dE1_rho
    H23 = 0 
    H33 = 0
    H34 = 0
    H35 = k
  } else if (poly_spec==3){
    H11 = 0
    H12 = 0
    H13 = 0 
    H14 = 0
    H15 = k
    H23 = 0  
    H33 = 0
    H34 = 0
    H35 = k  
  }
  
  ## line 4
  A4 = -t(weights) %*% t(x)
  ### element 4
  ### element 5
  H44 = A4 %*% dE1_beta
  H45 = A4 %*% dE1_rho
  
  ## line 5
  ### element 5
  H55 = -t(z) %*% dE1_rho
  
  Hess =  rbind( c(H11, H12, H13, H14, H15),
                 c(H12, H22, H23, H24, H25),
                 c(H13, H23, H33, H34, H35),
                 c(H14, H24, H34, H44, H45),
                 cbind(t(H15), t(H25), t(H35), t(H45), H55)) 
  #browser()
  return(Hess)   
}

# Model: ADL-MIDAS (FADL-MIDAS) regression 
# Estimation: MIDAS-NLS based on analytic derivatives (gradient and Hessian)
# Polynomial: Exponential Almon parametrized by theta1 and theta 2.
# Paper, analytic derivatives: based on  
# - Kostrov A. and Tetereva A. (2018) Forecasting realized correlations: a MIDAS approach. Working paper. 
# - https://www.researchgate.net/publication/331498346_Forecasting_realized_correlations_a_MIDAS_approach
# Optimizer: "nlminb" with constraints. 
optim_ardl_expalmon <- function(y, z, x, nbtrials = 100){
  n <- length(y)
  k <- dim(z)[2]
  if (is.null(k)){
    z <- matrix(z, nrow=length(z))
    k <- 1
  }
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)
  
  # store data into a list:
  dataList <- list (y=y, z=z, x=x)
  # get parameter constraints:
  constr <- get_constr_expalmon(k)
  
  #-------------------- main optimization --------------------#
  invisible(capture.output(opt <- multiStartoptim(objectivefn = estimate_ardl_expalmon,
                                        gradient = gradient_ardl_expalmon,
                                        data = dataList,
                                        hessian = hessian_ardl_expalmon,
                                        lower = constr$lw_b, upper = constr$up_b, 
                                        method = "nlminb",
                                        nbtrials = nbtrials,
                                        typerunif = "runifbase",
                                        control = list(eval.max=1000, iter.max=1000,
                                                       rel.tol=1e-12, x.tol=1.5e-10))))
  #-------------------- xxxxxxxxxxxxxxxx --------------------#
  # back-out parameters:
  k1 <- opt$par[1]
  k2 <- opt$par[2]
  beta <- opt$par[3]
  c <- opt$par[4]
  rho <- opt$par[5:length(opt$par)]

  
  
  # sort the output:
  other_x <- NULL
  for (i in 1:k)
    other_x <- c(other_x, paste0("other-",i))
  
  int <- "(Intercept)"
  coef <- matrix(c(c, rho, beta, k1, k2), nrow = 1, ncol = k + 4)
  rownames(coef) <- ""
  colnames(coef) <- c(int,other_x,"beta","k1","k2")
  class(coef) <- "midas.expalmon"
  return(coef)
}

get_constr_expalmon <- function(k){
  # define upper and lower constraints
  up_beta <- 20
  up_rho <- 20   
  up_theta1 <- 2   
  up_theta2 <- 0   
  
  low_beta <- -10
  low_rho <- -20     
  low_theta1 <-  0   
  low_theta2 <- -2 
  
  up_b <- c(up_theta1 , up_theta2, up_beta,  rep(up_rho, k+1))
  lw_b <- c(low_theta1, low_theta2, low_beta, rep(low_rho, k+1))  
  
  return(list(up_b = up_b, lw_b = lw_b)) 
}

estimate_ardl_expalmon <- function(args, dataList) {
  y <- dataList$y
  z <- dataList$z
  x <- dataList$x
  
  p <- dim(x)[2]
  xi <- matrix(1:p, p,1)
  xi_sq <- xi^2  
  
  theta1 <- args[1]  
  theta2 <- args[2]
  beta <- args[3]
  rho <- args[4:length(args)] 
  
  weights <- exp(theta1 * xi + theta2 * xi_sq)
  #browser()
  mse <- t(y  - z %*% rho -  beta * x %*% weights)%*%(y  - z %*% rho -  beta * x %*% weights)  
  
  return(as.numeric(mse))
  
}

gradient_ardl_expalmon <- function(args, dataList) {
  y <- dataList$y
  z <- dataList$z
  x <- dataList$x
  
  p <- dim(x)[2]
  xi <- matrix(1:p, p,1)
  xi_sq <- xi^2  
  
  theta1 <- args[1]  
  theta2 <- args[2]
  beta <- args[3]
  rho <- args[4:length(args)] 
  
  weights <- exp(theta1 * xi + theta2 * xi_sq) 
  
  
  nabla_G <- 2*(y- beta * x %*%  weights - z%*%rho )
  
  T_star_w_C <-  -beta * t(x)
  T_star_beta_C <-  -t(weights) %*% t(x)
  T_star_rho_C <-  - t(z)
  
  
  T_star_theta1_W <- t(weights * xi)
  T_star_theta2_W <- t(weights * xi_sq)
  
  
  part1 <- T_star_theta1_W %*% T_star_w_C %*% nabla_G   # theta1
  part2 <- T_star_theta2_W %*% T_star_w_C %*% nabla_G   # theta2
  part3 <- T_star_beta_C  %*% nabla_G   # beta
  part4 <- T_star_rho_C  %*% nabla_G    # rho
  grad <-  c(part1, part2, part3, part4)
  
  return(grad)  
}

hessian_ardl_expalmon <- function(args, dataList) {
  y <- dataList$y
  z <- dataList$z
  x <- dataList$x
  
  p <- dim(x)[2]
  xi <- matrix(1:p, p,1)
  xi_sq <- xi^2  
  
  theta1 <- args[1]  
  theta2 <- args[2]
  beta <- args[3]
  rho <- args[4:length(args)] 
  
  weights <- exp(theta1 * xi + theta2 * xi_sq) 
  
  ## line 1
  A1_1 <- t(weights * xi)
  A1_2 <- -beta * t(x)
  A1 <-  A1_1 %*% A1_2 
  B1 <- 2*(y - beta * x %*%  weights - z%*%rho )
  ### element 1
  dA1_theta1 <- t(weights * xi_sq) %*%   A1_2
  dB1_theta1 <-  - 2* beta * x %*%  (weights * xi)
  ### element 2
  dA1_theta2 <- t(weights * xi^3) %*%   A1_2
  dB1_theta2 <-  - 2* beta * x %*%  (weights * xi_sq)
  ### element 3
  dA1_beta <- t(weights * xi)%*%t(-x)
  dB1_beta <- -2*x%*%weights
  ### element 4
  dB1_rho <- -2*z
  H11 <- dA1_theta1 %*% B1 + A1 %*% dB1_theta1 
  H12 <- dA1_theta2 %*% B1 + A1 %*% dB1_theta2 
  H13 <- dA1_beta %*% B1 + A1 %*% dB1_beta 
  H14 <- A1 %*% dB1_rho
  
  ## line 2
  A2_1 <- t(weights * xi_sq)
  A2 <- A2_1 %*% A1_2 
  ### element 2
  dA2_theta2 <- t(weights * xi^4) %*%   A1_2
  ### element 3
  dA2_beta <- A2_1 %*% t(-x)
  ### element 4
  H22 <- dA2_theta2 %*% B1 +  A2 %*% dB1_theta2 
  H23 <- dA2_beta %*% B1 + A2 %*% dB1_beta 
  H24 <- A2 %*% dB1_rho
  
  ## line 3
  A3 <-  -t(weights) %*% t(x)
  ### element 3
  ### element 4
  ### element 5
  H33 <- A3  %*% dB1_beta
  H34 <- A3  %*% dB1_rho 
  
  
  ## line 4
  H44 <- t(-z) %*% dB1_rho 
  
  Hess <- rbind( c(H11, H12, H13, H14),
                 c(H12, H22, H23, H24),
                 c(H13, H23, H33, H34),
                 cbind(t(H14), t(H24), t(H34), H44)) 
  #browser()
  return(Hess)  
}

# Model: MIDAS Logit regression 
# Estimation: MIDAS-NLS based on analytic derivatives (gradient and Hessian)
# Polynomial: Exponential Almon parametrized by theta1 and theta 2
# Paper, analytic derivatives: based on  
# - 	F. Audrino, A. Kostrov, and J.-P. Ortega (2019) Predicting U.S. Bank Failures with MIDAS Logit Models, JFQA, 54(6), 2575-2603. 
# - https://doi.org/10.1017/S0022109018001308
# Optimizer: nlm.
optim_logit_expalmon <- function(y, z, x){
  if (!all(unique(y)%in%c(0,1))){
    stop("response variable must be binary for MIDAS Logit model.")
  }
  n <- length(y)
  k <- dim(z)[2]
  if (is.null(k)){
    z <- matrix(z, nrow=length(z))
    k <- 1
  }
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)

  # initial values
  init = rep(0,3 + k+1) 
  #-------------------- main optimization --------------------#
  opt <- suppressWarnings(stats::nlm(f = estimateLogit_expalm_nlm,
             p = init,
             y=y, V=z, Z=x,
             hessian = TRUE,
             check.analyticals = TRUE))
  #-------------------- xxxxxxxxxxxxxxxx --------------------#
  # back-out parameters:
  k1 <- opt$estimate[1]
  k2 <- opt$estimate[2]
  beta <- opt$estimate[3]
  c <- opt$estimate[4]
  rho <- opt$estimate[5:length(opt$estimate)]

  
  # sort the output:
  fn <- NULL
  for (i in 1:k)
    fn <- c(fn, paste0("factor-",i))

  int <- "(Intercept)"
  coef <- matrix(c(c, rho, beta, k1, k2), nrow = 1, ncol = k + 4)
  rownames(coef) <- ""
  colnames(coef) <- c(int,fn,"beta","k1","k2")
  class(coef) <- "logit.expalmon"
  return(coef)
}

estimateLogit_expalm_nlm <- function(args, y, V, Z){
  
  iota = matrix(1, length(y) )
  p = dim(Z)[2]
  xi = matrix(1:p, p,1)
  xi_sq = xi^2  
  
  theta1 = args[1]  
  theta2 = args[2]
  beta = args[3]
  rho = args[4:length(args)] 
  
  weights = exp(theta1 * xi + theta2 * xi_sq)
  
  
  LogLM = -( t(y) %*% (V %*% rho  + beta * Z %*% weights) - 
               t(iota) %*% log(iota+exp(V %*% rho  + beta * Z %*% weights))  )
  
  
  out = as.numeric(LogLM) 
  if (is.na(out)==TRUE) {
    out = Inf
  }
  
  attr(out, 'gradient') <-  gradLogit_expalm(args, y, V, Z) 
  attr(out, 'hessian') <-   HessianLogit_expalm(args, y, V, Z)
  
  return(out)
}

gradLogit_expalm <- function(args, y, V, Z) {
  
  p = dim(Z)[2]
  xi = matrix(1:p, p,1)
  xi_sq = xi^2  
  
  theta1 = args[1]  
  theta2 = args[2]
  beta = args[3]
  rho = args[4:length(args)] 
  
  weights = exp(theta1 * xi + theta2 * xi_sq) 
  
  nabla_G = - ( y- exp(beta * Z %*%  weights + V%*%rho) *
                  1 / (1 + exp(beta * Z %*%  weights + V%*%rho)))
  
  T_star_w_C =  beta * t(Z)
  T_star_beta_C =  t(weights) %*% t(Z)
  T_star_rho_C =   t(V)
  
  T_star_theta1_W = t(weights * xi)
  T_star_theta2_W = t(weights * xi_sq)
  
  part1 = T_star_theta1_W %*% T_star_w_C %*% nabla_G   # theta1
  part2 = T_star_theta2_W %*% T_star_w_C %*% nabla_G   # theta2
  part3 = T_star_beta_C  %*% nabla_G   # beta
  part4 = T_star_rho_C  %*% nabla_G    # rho
  grad =  c(part1, part2, part3, part4)
  
  return(grad)  
}


HessianLogit_expalm <- function(args, y, V, Z) {
  
  p = dim(Z)[2]
  xi = matrix(1:p, p,1)
  xi_sq = xi^2  
  
  theta1 = args[1]  
  theta2 = args[2]
  beta = args[3]
  rho = args[4:length(args)] 
  
  weights = exp(theta1 * xi + theta2 * xi_sq) 
  lincomb =  beta * Z %*%  weights + V%*%rho 
  
  ## line 1
  A1_1 = t(weights * xi)
  A1_2 = beta * t(Z)
  A1 =  A1_1 %*% A1_2 
  B1 = - ( y- exp(lincomb) * 1 / (1 + exp(lincomb)))
  dB1_lincomb = exp(lincomb)*(1+exp(lincomb))^(-1) - exp(lincomb)^(2) *  (1+exp(lincomb))^(-2)
  
  ### element 1
  dA1_theta1 = t(weights * xi_sq) %*%   A1_2
  dB1_theta1 =  dB1_lincomb * (beta * Z) %*% (weights * xi)       
  ### element 2
  dA1_theta2 = t(weights * xi^3) %*%   A1_2
  dB1_theta2 = dB1_lincomb * (beta * Z) %*% (weights * xi_sq) 
  ### element 3
  dA1_beta = t(weights * xi)%*%t(Z)
  dB1_beta =   dB1_lincomb * (Z%*%weights)
  ### element 4
  dB1_rho =   diag( as.vector(dB1_lincomb)  ) %*% V
  H11 = dA1_theta1 %*% B1 + A1 %*% dB1_theta1 
  H12 = dA1_theta2 %*% B1 + A1 %*% dB1_theta2 
  H13 = dA1_beta %*% B1 + A1 %*% dB1_beta 
  H14 = A1 %*% dB1_rho
  
  ## line 2
  A2_1 = t(weights * xi_sq)
  A2 = A2_1 %*% A1_2 
  ### element 2
  dA2_theta2 = t(weights * xi^4) %*%   A1_2
  ### element 3
  dA2_beta = A2_1 %*% t(Z)
  ### element 4
  H22 = dA2_theta2 %*% B1 +  A2 %*% dB1_theta2 
  H23 = dA2_beta %*% B1 + A2 %*% dB1_beta 
  H24 = A2 %*% dB1_rho
  
  ## line 3
  A3 =  t(weights) %*% t(Z)
  ### element 3
  ### element 4
  ### element 5
  H33 = A3  %*% dB1_beta
  H34 = A3  %*% dB1_rho 
  
  ## line 4
  H44 = t(V) %*% dB1_rho 
  
  Hess = rbind( c(H11, H12, H13, H14),
                c(H12, H22, H23, H24),
                c(H13, H23, H33, H34),
                cbind(t(H14), t(H24), t(H34), H44)) 
  #browser()
  return(Hess)  
}


optim_ardl_legendre <- function(y, z, x, legendre_degree = 3){
  n <- length(y)
  k <- dim(z)[2]
  if (is.null(k)){
    z <- matrix(z, nrow=length(z))
    k <- 1
  }
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)
  p <- lb(degree = legendre_degree, jmax = dim(x)[2])
  x <- x%*%p
  bigx <- cbind(z, x)
  
  coef <- t(solve(t(bigx)%*%bigx)%*%t(bigx)%*%y)
  other_x <- NULL
  for (i in 1:k) 
    other_x <- c(other_x, paste0("other-",i))
  
  poly_n <- NULL
  for (i in 0:legendre_degree) 
    poly_n <- c(poly_n, paste0("Poly-",i))
  
  int <- "(Intercept)"
  rownames(coef) <- ""
  colnames(coef) <- c(int,other_x,poly_n)
  class(coef) <- "midas.legendre"
  return(coef)
}


optim_logit_legendre <- function(y, z, x, legendre_degree = 3){
  n <- length(y)
  k <- dim(z)[2]
  if (is.null(k)){
    z <- matrix(z, nrow=length(z))
    k <- 1
  }
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)
  p <- lb(degree = legendre_degree, jmax = dim(x)[2])
  x <- x%*%p
  bigx <- cbind(z, x)
  
  coef <- stats::glm(y~-1+bigx)$coef
  coef <- matrix(coef, ncol = length(coef))
  fn <- NULL
  for (i in 1:k)
    fn <- c(fn, paste0("factor-",i))
  
  poly_n <- NULL
  for (i in 0:legendre_degree) 
    poly_n <- c(poly_n, paste0("Poly-",i))
  
  int <- "(Intercept)"
  rownames(coef) <- ""
  colnames(coef) <- c(int,fn,poly_n)
  class(coef) <- "logit.legendre"
  return(coef)
}

multiStartoptim <-  function(start0 = NULL, objectivefn,  gradient = NULL, ..., hessian = NULL, 
                             lower = -Inf, upper = Inf, control = list(),
                             method = c("L-BFGS-B", "Nelder-Mead","nlminb"), nbtrials = NULL, 
                             typerunif = c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp"), localsearch = c("exhaustive", "median"), 
                             verb=FALSE, nbclusters = NULL)
{  
  if(missing(objectivefn)) stop("Objective function must be provided")
  
  nbpar <- length(lower)       
  if (nbpar != length(upper) || lower > upper) stop("parameter lower must be lower than parameter upper componentwise, and they must have the same length")
  
  if (missing(method) || !(method %in% c("L-BFGS-B", "Nelder-Mead", "nlminb")))  
  {
    warning("optimization method is missing or inadequate, default is set to nlminb")
    themethod <- "nlminb"
  }
  else  
  {
    themethod <- match.arg(method)
  }
  
  if (missing(nbclusters) || is.null(nbclusters) || nbclusters <= 1)
  {      
    if (missing(start0) && (missing(nbtrials) || nbtrials <= 0)) stop("When nbtrials is missing or equal to 0, starting parameters 'start0' should be provided")
    
    if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0)
    {        
      if (!missing(typerunif)) warning("unused argument typerunif")
      
      if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
      
      
      return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                             lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                             verb = verb, cl = NULL))          
    } 
    else
    {          
      if (!missing(start0)) warning("starting value is unused when nbtrials is provided")
      
      if (missing(typerunif) || !(typerunif %in% c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp")))
      {
        warning("typerunif is either missing or inadequate, and was set to default 'runifbase'")
        typerunif <- "runifbase"
      }
      else 
      {
        typerunif <- match.arg(typerunif)
      }
      
      if (!is.finite(lower) || !is.finite(upper))
        stop("Both the lower and upper bounds for parameters should be provided and finite")
      
      localsearch <- match.arg(localsearch)
      
      if (length(localsearch) == 0) localsearch <- "exhaustive"
      
      if (localsearch == "exhaustive")
      {
        start0 <- rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                        method=typerunif)
      }
      
      if (localsearch == "median")
      {                                    
        if (nbpar == 1)
        {  
          U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                      method=typerunif)
          fwU <- objectivefn(U, ...)
          start0 <- U[fwU < median(fwU)] 
          nbtrials <- length(start0)
        }
        else
        { 
          U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                      method=typerunif)
          fwU <- apply(U, 1, function(x) objectivefn(x, ...))
          start0 <- U[fwU < median(fwU), ]
          nbtrials <- dim(start0)[1]
        }    
      }
      
      if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
      return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                             lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                             verb = verb, cl = NULL))                 
    } 
  } 
  else
  {
    if (nbclusters < 0 || floor(nbclusters) != nbclusters || is.array(nbclusters) || !is.numeric(nbclusters)) stop("nbclusters must be a positive integer")
    
    if(missing(lower) || missing(upper)) stop("lower and upper are missing and must be provided")
    
    if (!is.finite(lower) || !is.finite(upper)) stop("Both the lower and upper bounds for parameters should be provided and finite")
    
    if(verb) warning("argument verb is unused in parallel computation")
    
    if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0) 
    {
      warning("nbtrials is not provided (or negative), default is set to 50")
      nbtrials <- 50
    }
    
    if (missing(typerunif) || !(typerunif %in% c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp")))
    {
      warning("typerunif is either missing or inadequate, and was set to default 'runifbase'")
      typerunif <- "runifbase"
    } else 
    {
      typerunif <- match.arg(typerunif)
    }
    
    localsearch <- match.arg(localsearch)
    
    if (length(localsearch) == 0) localsearch <- "exhaustive"
    
    if (localsearch == "exhaustive")
    {
      start0 <- rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                      method=typerunif)
    }
    
    if (localsearch == "median")
    {                                    
      if (nbpar == 1)
      {  
        U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                    method=typerunif)
        fwU <- objectivefn(U, ...)
        start0 <- U[fwU < median(fwU)] 
        nbtrials <- length(start0)
      }
      else
      { 
        U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                    method=typerunif)
        fwU <- apply(U, 1, function(x) objectivefn(x, ...))
        start0 <- U[fwU < median(fwU), ]
        nbtrials <- dim(start0)[1]
      }    
    }     
    
    if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
    
    packageStartupMessage("Processing...", appendLF = FALSE)  
    cat("\n")    
    return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                           lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                           verb = verb, cl = nbclusters))
  }
  
}

rUnif <- function(nbtrials, nbpar, lower, upper, 
                  method=c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp"))
{
  callrunifmat <- function(nbtrials, nbpar, lower, upper)
  { 
    if (nbpar == 1) 
    {
      return(lower + (upper - lower) * randtoolbox::SFMT(nbtrials))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::SFMT(nbtrials,nbpar))
    }
  }
  
  callrunifmatantith <- function(nbtrials, nbpar, lower, upper)
  { 
    
    nbtrials_2 <- 0.5*nbtrials
    if (nbpar == 1) 
    {
      Uni <- randtoolbox::SFMT(nbtrials_2)
      return(lower + (upper - lower) * c(Uni, 1-Uni))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
      etendue <- upper - lower
      Uni <- randtoolbox::SFMT(nbtrials_2,nbpar)
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol = nbpar, byrow = T)*rbind(Uni, 1-Uni))
    }
  }
  
  callxlodisp <- function(nbtrials, lower, upper)
  {
    nbpar <- length(lower)
    
    if (nbpar == 1) 
    { 
      ylodisp <- log(2*seq_len(nbtrials)[-1] - 3)/(log(2)) 
      xlodisp <- c(1, ylodisp - floor(ylodisp))     
      return(lower + (upper - lower) * xlodisp)
    }
    else
    {
      ylodisp <- log(2*seq_len(nbtrials)[-1] - 3)/(log(2)) 
      xlodisp <- c(1, ylodisp - floor(ylodisp))
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)*
               sapply(seq_len(nbpar), function (x) sample(xlodisp, nbtrials, replace = FALSE)))
    }
  }
  
  if (method == "runifbase")
  {
    return (callrunifmat(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper))
  }
  
  if (method == "runifantithetics")
  {
    if(floor(nbtrials) != nbtrials) 
    {
      stop("To achieve this, nbtrials should be an even number")
    }    
    return(callrunifmatantith(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper))
  }
  
  if (method == "sobol")
  {
    if (nbpar == 1) 
    {
      return(lower + (upper - lower) * randtoolbox::sobol(n = nbtrials, dim = nbpar, scrambling = 2))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::sobol(n = nbtrials, dim = nbpar, scrambling = 2))
    }
  }
  
  if (method == "torus")
  {
    if (nbpar == 1) 
    {
      return(lower + (upper - lower) * randtoolbox::torus(n = nbtrials, dim = nbpar))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::torus(n = nbtrials, dim = nbpar))
    }
  }
  
  if (method == "niederreiterlodisp")
  {
    return(callxlodisp(nbtrials = nbtrials, lower = lower, upper = upper))
  }
}


call_localoptim <- function(start0, objectivefn,  gradient = NULL, ..., hessian = NULL, 
                            lower, upper, control = list(), themethod, nbtrials = NULL, 
                            verb=FALSE, cl = NULL)
{ 
  nbpar <- length(lower)
  
  if (missing(cl) || is.null(cl))
  {
    if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0)
    {            
      if (themethod == "L-BFGS-B")
      {
        if (missing(hessian)) hessian <- FALSE
        return(optim(par = start0, fn = objectivefn, gr = gradient,  ..., 
                     method = themethod, lower = lower, upper = upper, 
                     control = control, hessian = hessian))            
      } 
      
      if (themethod == "Nelder-Mead")
      {
        if (missing(hessian)) hessian <- FALSE
        return(optim(par = start0, fn = objectivefn, gr = gradient,  ..., 
                     method = themethod, control = control, hessian = hessian))
      }
      
      if (themethod == "nlminb") 
      { 
        
        return(nlminb(start = start0, objective = objectivefn, gradient = gradient, 
                      hessian = hessian, ..., scale = 1, control = control, 
                      lower = lower, upper = upper))                    
      }
    } 
    else 
    {        
      ObjectifNSScour <- 1000000
      ObjectifNSSi <- 0
      res_cour <- 0 
      nbpar <- length(lower)
      if (length(upper) != nbpar || lower > upper) stop("lower and upper are incorrect")
      
      if (verb == TRUE)
      {            
        verbtrace <- list(par = matrix(NA, nrow=nbtrials, ncol=nbpar), 
                          objective = rep.int(NA, nbtrials),
                          iteration = rep.int(NA, nbtrials))
      }
      
      howFar <- txtProgressBar(min=0,max=nbtrials,style=3)
      
      for(i in seq_len(nbtrials))
      {          
        if(nbpar == 1)
        {
          start_0 <- start0[i]
        } else              
        {
          start_0 <- start0[i,]
        }
        
        if (themethod == "L-BFGS-B")
        { 
          if (missing(hessian)) hessian <- FALSE
          res_cour <- suppressWarnings(try(optim(par = start_0, fn = objectivefn, gr = gradient,  ..., 
                                                 method = themethod, lower = lower, upper = upper, 
                                                 control = control, hessian = hessian), silent = TRUE))
          
          try(ObjectifNSSi <- res_cour$value, silent = TRUE)
          
        } 
        
        if (themethod == "Nelder-Mead")
        { 
          if (missing(hessian)) hessian <- FALSE
          res_cour <- suppressWarnings(try(optim(par = start_0, fn = objectivefn, gr = gradient,  ..., 
                                                 method = themethod, 
                                                 control = control, hessian = hessian), silent = TRUE))
          
          
          ObjectifNSSi <- res_cour$value
          
        }
        
        if (themethod == "nlminb") 
        { 
          res_cour <- suppressWarnings(try(nlminb(start = start_0, objective = objectivefn, gradient = gradient, hessian = hessian,
                                                  ..., scale = 1, control = control, lower = lower, upper = upper), silent = TRUE))
          
          
          ObjectifNSSi <- res_cour$objective
          
        }
        
        if(ObjectifNSSi < ObjectifNSScour) 
        {
          
          res <- res_cour
          ObjectifNSScour <- ObjectifNSSi
          if (verb == TRUE)
          {
            verbtrace$par[i, ] <- res_cour$par
            verbtrace$objective[i] <- ObjectifNSSi
            verbtrace$iteration[i] <- i
          }
        }
        
        
        setTxtProgressBar(howFar, value=i)
      }    
      close(howFar)
      
      if (verb == TRUE)
      { 
        booltrace <- !is.na(verbtrace$iteration)
        if (nbpar == 1)
        {startingparams <- start0[booltrace]}
        else startingparams <- start0[booltrace,]
        foundparams <- verbtrace$par[booltrace, ]
        objective_val <- verbtrace$objective[booltrace]
        iteration_no <- verbtrace$iteration[booltrace]
        return(list(res=res, 
                    iteration_no = iteration_no,
                    startingparams_sequence = startingparams,
                    foundparams_sequence = foundparams, 
                    objective_val_sequence = objective_val))  
      } else return(res)        
    }
  } 
  else 
  {
    
    force(objectivefn)     
    if (!missing(gradient)) force(gradient)    
    
    if (length(upper) != nbpar || lower > upper) stop("lower and upper are incorrect")
    
    cl <- makeCluster(cl, type="SOCK")    
    .env <- environment()
    ldots <- list(...)  
    clusterExport(cl, "ldots", envir = .env)
    
    if (themethod == "L-BFGS-B")
    {       
      if (missing(hessian)) hessian <- FALSE
      if (nbpar <= 1)
      {
        objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                      method = themethod,  
                                                                                      lower = lower, upper = upper,
                                                                                      control = control, hessian = hessian), silent=TRUE)))
        on.exit(stopCluster(cl))           
        packageStartupMessage("             Done !")        
        
        return(objectivemin[ , which.min(objectivemin[2,])])
      } 
      else 
      {
        objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                      method = themethod,  
                                                                                      lower = lower, upper = upper,
                                                                                      control = control, hessian = hessian), silent=TRUE)))
        on.exit(stopCluster(cl))
        packageStartupMessage("             Done !")        
        
        return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "value")))]])
        
      }
    }
    
    if (themethod == "Nelder-Mead")
    {   
      
      if (missing(hessian)) hessian <- FALSE
      if (nbpar <= 1)
      {
        objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                      method = themethod,
                                                                                      control = control, hessian = hessian), silent=TRUE)))
        on.exit(stopCluster(cl))
        packageStartupMessage("             Done !")        
        
        return(objectivemin[ , which.min(objectivemin[2,])])
      } 
      else 
      {
        objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                      method = themethod,
                                                                                      control = control, hessian = hessian), silent=TRUE)))
        on.exit(stopCluster(cl))
        packageStartupMessage("             Done !")        
        
        return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "value")))]])       
      }
    }
    
    if (themethod == "nlminb")
    {      
      if (!missing(hessian)) force(hessian)
      
      if (nbpar <= 1)
      {
        objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(nlminb(x, objective = objectivefn,
                                                                                       gradient = gradient, hessian = hessian, ...
                                                                                       , scale = 1, control = control,
                                                                                       lower = lower, upper = upper), silent=TRUE)))
        on.exit(stopCluster(cl))  
        packageStartupMessage("             Done !")        
        
        return(objectivemin[ , which.min(objectivemin[2,])])
      } 
      else
      {
        objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(nlminb(x, objective = objectivefn,
                                                                                       gradient = gradient, hessian = hessian, ...
                                                                                       , scale = 1, control = control,
                                                                                       lower = lower, upper = upper), silent=TRUE)))
        on.exit(stopCluster(cl))
        packageStartupMessage("             Done !")        
        
        return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "objective")))]])
        
      }    
    }
    
  }
}