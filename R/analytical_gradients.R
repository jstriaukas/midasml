# Model: ADL-MIDAS (FADL-MIDAS) regression 
# Estimation: MIDAS-NLS based on analytic derivatives (gradient and Hessian)
# Polynomial: Beta, 4 parametrizations, refer to "poly_spec"
# Paper, analytic derivatives:  
# - Kostrov (2020) Estimating MIDAS regressions via MIDAS-NLS with revised optimization. Working paper. 
# Optimizer: "nlminb", with constraints. 
optim_ardl_beta <- function(y, z, x, poly_spec, num.coef){
  n <- length(y)
  k <- dim(z)[2]
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)
  # store data into a list:
  dataList <- list (y=y, z=z, x=x, poly_spec=poly_spec)
  # get parameter constraints:
  constr <- get_constr_beta(poly_spec, k)
  
  #-------------------- main optimization --------------------#
  opt <- mcGlobaloptim::multiStartoptim(objectivefn = estimate_ardl_beta, 
                                        gradient = gradient_ardl_beta,
                                        data = dataList,
                                        hessian = hessian_ardl_beta, 
                                        lower = constr$lw_b, upper = constr$up_b, 
                                        method = "nlminb", 
                                        nbtrials = num.coef,
                                        typerunif = "runifbase",
                                        control = list(eval.max=1000, iter.max=1000,
                                                       rel.tol=1e-12, x.tol=1.5e-10))
  #-------------------- xxxxxxxxxxxxxxxx --------------------#
  # back-out parameters:
  k1 <- opt$par[1]
  k2 <- opt$par[2]
  k3 <- opt$par[3]
  beta <- opt$par[4]
  c <- opt$par[5]
  rho <- opt$par[6:length(opt$par)]
  
  # sort the output:
  ar_lags <- NULL
  for (i in 1:k) 
    ar_lags <- c(ar_lags, paste0("AR-",i))
  
  if (poly_spec == 0) {
    coef <- matrix(c(c, rho, beta, k1, k2, k3), nrow = 1, ncol = k + 5)
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,"beta","k1","k2", "k3")
  } else if (poly_spec == 1) {
    coef <- matrix(c(c, rho, beta, k2, k3), nrow = 1, ncol = k + 4)
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,"beta","k2", "k3")
  } else if (poly_spec == 2) {
    coef <- matrix(c(c, rho, beta, k1, k2), nrow = 1, ncol = k + 4)
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,"beta","k1", "k2")
  }
  return(coef)
}



get_constr_beta <- function(poly_spec, k){
  # define upper and lower constraints
  up_theta1 <- 2
  up_theta2 <- 20
  up_theta3 <- 0.7
  up_beta <- 5
  up_rho <- 7
  
  low_theta1 <- 0.1
  low_theta2 <- 1.2
  low_theta3 <- 0
  low_beta <- -5
  low_rho <- -7 
  
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
optim_ardl_expalmon <- function(y, z, x, num.coef){
  n <- length(y)
  k <- dim(z)[2]
  # append a vector of ones:
  iota <- matrix(1, nrow = n, ncol = 1)
  z <- cbind(iota, z)
  # store data into a list:
  dataList <- list (y=y, z=z, x=x)
  # get parameter constraints:
  constr <- get_constr_expalmon(k)
  
  #-------------------- main optimization --------------------#
  opt <- mcGlobaloptim::multiStartoptim(objectivefn = estimate_ardl_expalmon,
                                        gradient = gradient_ardl_expalmon,
                                        data = dataList,
                                        hessian = hessian_ardl_expalmon,
                                        lower = constr$lw_b, upper = constr$up_b, 
                                        method = "nlminb",
                                        nbtrials = num.coef,
                                        typerunif = "runifbase",
                                        control = list(eval.max=1000, iter.max=1000,
                                                       rel.tol=1e-12, x.tol=1.5e-10))
  #-------------------- xxxxxxxxxxxxxxxx --------------------#
  # back-out parameters:
  k1 <- opt$par[1]
  k2 <- opt$par[2]
  beta <- opt$par[3]
  c <- opt$par[4]
  rho <- opt$par[5:length(opt$par)]
  
  # sort the output:
  ar_lags <- NULL
  for (i in 1:k)
    ar_lags <- c(ar_lags, paste0("AR-",i))
  
  coef <- matrix(c(c, rho, beta, k1, k2), nrow = 1, ncol = k + 4)
  rownames(coef) <- ""
  colnames(coef) <- c("(Intercept)",ar_lags,"beta","k1","k2")
  
  return(coef)
}

get_constr_expalmon <- function(k){
  # define upper and lower constraints
  up_beta <- 3
  up_rho <- 10   
  up_theta1 <- 0.2   
  up_theta2 <- 0   
  
  low_beta <- -3
  low_rho <- -10     
  low_theta1 <-  0   
  low_theta2 <- -0.3 
  
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