// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <iostream>

#include <math.h>
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


//' Computes fast OLS solution
//' 
//' @param Y (vec) vector of response.
//' @param X (mat) matrix of covariates.
//' @param intercept (double) 1 - includes intercept, else not.
//' @export fastols
//' @keywords internal
// [[Rcpp::export]]
arma::vec fastols(const arma::vec & Y, const arma::mat & X, double intercept) {
  arma::mat X1;
  int n = X.n_rows;
  arma::mat ones(n, 1, fill::ones);
  if (intercept==1){
    X1 = join_rows(ones, X );
  } else {
    X1 = X;
  }
  arma::vec coef = arma::solve(X1, Y); 
  return coef;
}

//' Computes fast ALS solution
//' 
//' @param Y (vec) vector of response.
//' @param X (mat) matrix of covariates.
//' @param intercept (double) 1 - includes intercept, else not.
//' @param tau (double) quantile level.
//' @param maxIter maximum number of iterations.
//' @param thresh threshold value for stopping criterion.
//' @export fastals
//' @keywords internal
// [[Rcpp::export]]
arma::vec fastals(const arma::vec & Y, const arma::mat & X, double intercept, double tau, double maxIter, double thresh) {
  //------------------------------- fastals -----------------------------------------------------//
  // Purpose: compute fast Asymmetric Least Squares by using iteratively reweighted least squares.
  //---------------------------------------------------------------------------------------------//
  //add intercept if needed:
  arma::mat X1;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec beta(p);
  beta.fill(datum::nan);
  if (n<p){
    ::Rf_warning("n smaller than p. use penalized regression to compute the solution."); 
  } else {
    arma::mat ones(n, 1, fill::ones);
    if (intercept==1){
      X1 = join_rows(ones, X );
    } else {
      X1 = X;
    }
    // initialize variables used in iterations:
    arma::mat W;
    arma::vec xtwy;
    arma::mat xtwx, inv_xtwx;
    
    // compute OLS solution to initialize:
    beta = arma::solve(X1, Y); 
    arma::vec beta_old = beta;
    
    // iterated weighted least squares:
    double iter = 0;
    while (iter < maxIter){
      
      arma::vec xbeta = X1*beta;
      arma::vec indicator(n, fill::zeros);
      for (int j = 0; j < n; j++){
        if (Y(j)<=xbeta(j)){
          indicator(j) = 1;
        }
      }
      arma::vec w = tau - indicator;
      arma::vec aw = abs(w);
      
      // compute weights
      W = diagmat(aw);
      
      // compute weighted XtX 
      xtwx = X1.t()*W*X1;
      // compute inverse of XtWX
      inv_xtwx = arma::inv(xtwx);
      // compute weighted XtY
      xtwy = X1.t()*W*Y;
      // compute beta
      beta = inv_xtwx*xtwy;
      
      // check convergence
      double eps = norm(beta-beta_old,2);
      if (eps <=thresh){
        // converged
        iter =  maxIter;
      } else {
        // update previous beta and counter
        beta_old = beta;
        iter++;
      }
    }
  }
  return beta;
}

arma::mat boundc(colvec x, colvec dx){
  int n = x.n_rows;
  arma::vec b = 1e20 + 0 * x;
  arma::uvec f = find(dx < 0);
  b(f) = -x(f)/dx(f);
  arma::mat bout(b);
  bout.reshape(n, 1);
  return(bout);
}

//' Computes fast rq solution
//' 
//' @param Y (vec) vector of response.
//' @param X (mat) matrix of covariates.
//' @param intercept (double) 1 - includes intercept, else not.
//' @param tau (double) quantile level.
//' @export fastrq
//' @keywords internal
// [[Rcpp::export]]
arma::vec fastrq(const arma::vec & Y, const arma::mat & X, double intercept, double tau) {
  //------------------------------- fastals -----------------------------------------------------//
  // Purpose: compute fast Quantile Regression by using interior point (Frisch-Newton) algorithm.
  //---------------------------------------------------------------------------------------------//
  //add intercept if needed:
  arma::mat X1;
  double t = X.n_rows;
  int p = X.n_cols;
  arma::vec a, b, x;
  arma::mat A, c, s, Q, AQ;
  double beta_solver, small, max_it, m=0, n=0, gap=0;
  arma::colvec q, rhs, dx, ds, dxdz, dsdw, xinv, sinv,xi;
  arma::mat fx, fs, fw, fz;
  arma::mat fp, fd;
  arma::rowvec r, y, z, w, dy, dz, dw;
  double mu=0, g=0;
  if (t < p){
    ::Rf_warning("n smaller than p. use penalized regression to compute the solution."); 
  } else {
    arma::mat ones(t, 1, fill::ones);
    if (intercept==1){
      X1 = join_rows(ones, X );
    } else {
      X1 = X;
    }
    
    arma::vec u(t, fill::ones);
    a = (1-tau)*u;
    
    // setup variables for LP:
    A = X1.t();
    c = -Y.t();
    b =  X1.t()*a;
    x = a;
    
    // some constants
    beta_solver = 0.9995;
    small = 1e-5;
    max_it = 50;
    m = A.n_rows;
    n = A.n_cols;
    
    // Generate inital feasible point
    s = u - x;
    y = solve(A.t(),c.t()).t();
    r = c - y * A;
    arma::uvec zero_idx = find(r == 0); 
    arma::uvec nzero_idx = find(r != 0); 
    arma::rowvec rtmp = r;
    rtmp.elem(zero_idx).fill(1.0);  
    rtmp.elem(nzero_idx).fill(0.0); 
    r = r + 0.001 * rtmp;   
    arma::uvec nonpos = find(r <= 0 ); 
    arma::uvec pos = find(r > 0); 
    arma::rowvec rtmp2 = r;
    rtmp2.elem(nonpos).fill(0.0);  
    rtmp2.elem(pos).fill(1.0); 
    z = r % rtmp2;
    w = z - r;
    arma::mat tmp = c * x - y * b + w * u;
    gap = tmp(0,0);
    
    // Start iterations
    int it = 0;
    while ((gap > small) && (it < max_it)){
      it = it + 1;
      //Compute affine step
      q = z.t()/x + w.t()/s;
      q =  1/q;
      r = z - w;
      Q = diagmat(sqrt(q));
      AQ = A * Q;         
      rhs = Q * r.t();        
      dy = solve(AQ.t(),rhs).t();
      dx = q % (dy * A - r).t();
      ds = -dx;
      dz = -z % (1 + dx / x).t();
      dw = -w % (1 + ds / s).t();
      // Compute maximum allowable step lengths
      fx = boundc(x, dx);
      fs = boundc(s, ds);
      fw = boundc(w.t(), dw.t());
      fz = boundc(z.t(), dz.t());
      fp = min(fx, fs);
      fd = min(fw, fz);
      fp = min(min(beta_solver * fp), 1);
      fd = min(min(beta_solver * fd), 1);
      
      // If full step is feasible, take it. Otherwise modify it
      double check;
      if (fp(0,0) > fd(0,0)){
        check = fd(0,0);
      } else {
        check = fp(0,0);
      }
      if (check < 1){
        // Update mu
        arma::mat tmp2 =  z * x + w * s;
        mu = tmp2(0,0);
        tmp2  = (z + fd(0,0) * dz)*(x + fp(0,0) * dx) + (w + fd(0,0) * dw) * (s + fp(0,0) * ds); 
        g = tmp2(0,0);
        double gmu = g/mu;
        mu = mu*pow(gmu,3)/( 2* n);
        
        // Compute modified step
        dxdz = dx % dz.t();
        dsdw = ds % dw.t();
        xinv = 1 / x;
        sinv = 1 / s;
        xi =  mu*(xinv - sinv);  
        rhs = rhs + Q*(dxdz - dsdw - xi);
        dy = solve(AQ.t(),rhs).t();
        dx = q % (A.t() * dy.t() + xi - r.t() -dxdz + dsdw);
        ds = -dx;
        dz = mu * xinv.t() - z - xinv.t() % z % dx.t() - dxdz.t();
        dw = mu * sinv.t() - w - sinv.t() % w % ds.t() - dsdw.t();
        
        // Compute maximum allowable step lengths
        fx = boundc(x, dx);
        fs = boundc(s, ds);
        fw = boundc(w.t(), dw.t());
        fz = boundc(z.t(), dz.t());
        fp = min(fx, fs);
        fd = min(fw, fz);
        fp = min(min(beta_solver * fp), 1);
        fd = min(min(beta_solver * fd), 1);
      }
      // Take the step
      x = x + fp(0,0) * dx;
      s = s + fp(0,0) * ds;
      y = y + fd(0,0) * dy;
      w = w + fd(0,0) * dw;
      z = z + fd(0,0) * dz;
      tmp = c * x - y * b + w * u;
      gap = tmp(0,0);
    }
  }
  arma::vec beta = -y.t();
  return beta;
}

//' Computes fast DL-MIDAS profiling solution
//' 
//' @param Y (vec) vector of response.
//' @param X (mat) matrix of covariates.
//' @param intercept (double) 1 - includes intercept, else not.
//' @param tau (double) quantile level.
//' @param which_loss loss function choice. 1 - mse, 2 - als, 3 - rq.
//' @param num_evals number of evalution points of MIDAS parameter.
//' @export midas_pr
//' @keywords internal
// [[Rcpp::export]]
arma::vec midas_pr(const arma::vec & Y, const arma::mat & X, double intercept, double tau, double which_loss, double num_evals) {
  int n = X.n_rows, dlag = X.n_cols;
  arma::vec u = linspace<vec>(0, 1, dlag);
  arma::vec param = linspace<vec>(1, 100, num_evals);
  arma::vec w;
  arma::vec Xw;
  arma::vec b;
  arma::vec res;
  arma::vec fit(num_evals);
  fit.fill(datum::nan);
  arma::mat ones(n, 1, fill::ones);
  arma::mat Xw1;
  arma::mat tmp;
  for (int i=0; i<num_evals; i++){
    w = pow(1-u,param(i)-1);
    w = w/sum(w);
    Xw = X*w;
    if (which_loss==1){
      b = fastols(Y, Xw, intercept);
    }
    if (which_loss==2){
      b = fastals(Y, Xw, intercept, tau, 100, 1e-3);
    }
    if (which_loss==3){
      b = fastrq(Y, Xw, intercept, tau);
    }
    if (intercept==1){
      Xw1 = join_rows(ones, Xw);
    } else {
      Xw1 = Xw;
    }
    res = Y-Xw1*b;
    if (which_loss==1){
      tmp = sum(pow(res,2))/n;
    }
    if (which_loss==2){
      arma::vec ind(n, fill::zeros);
      arma::vec ind2(n, fill::ones);
      arma::uvec idx = find(res < 0);
      ind(idx) = ind2(idx);
      tmp = sum(pow(res,2)%(tau - ind))/n;
    }
    if (which_loss==3){
      arma::vec ind(n, fill::zeros);
      arma::vec ind2(n, fill::ones);
      arma::uvec idx = find(res < 0);
      ind(idx) = ind2(idx);
      tmp = sum(res%(tau - ind))/n;
    }
    fit(i) = tmp(0,0);
  }
  arma::uvec idx2 = find(fit == min(fit));
  arma::vec par = param(idx2);
  w = pow(1-u,par(0)-1);
  w = w/sum(w);
  Xw = X*w;
  if (which_loss==1){
    b = fastols(Y, Xw, intercept);
  }
  if (which_loss==2){
    b = fastals(Y, Xw, intercept, tau, 1000, 1e-7);
  }
  if (which_loss==3){
    b = fastrq(Y, Xw, intercept, tau);
  }
  arma::vec beta = join_cols(b,par);
  
  return(beta);
}

//' Computes fast ARDL-MIDAS profiling solution
//' 
//' @param Y (vec) vector of response.
//' @param YLAG (mat) matrix of lagged values of Y.
//' @param X (mat) matrix of covariates.
//' @param intercept (double) 1 - includes intercept, else not.
//' @param tau (double) quantile level.
//' @param which_loss loss function choice. 1 - mse, 2 - als, 3 - rq.
//' @param num_evals number of evalution points of MIDAS parameter.
//' @export midasar_pr
//' @keywords internal
// [[Rcpp::export]]
arma::vec midasar_pr(const arma::vec & Y, const arma::mat & YLAG, const arma::mat & X, double intercept, double tau, double which_loss, double num_evals) {
  int n = X.n_rows, dlag = X.n_cols;
  arma::vec u = linspace<vec>(0, 1, dlag);
  arma::vec param = linspace<vec>(1, 100, num_evals);
  arma::vec w;
  arma::vec Xw;
  arma::vec b;
  arma::vec res;
  arma::vec fit(num_evals);
  fit.fill(datum::nan);
  arma::mat ones(n, 1, fill::ones);
  arma::mat Xw1;
  arma::mat Xw_ar;
  arma::mat tmp;
  for (int i=0; i<num_evals; i++){
    w = pow(1-u,param(i)-1);
    w = w/sum(w);
    Xw = X*w;
    Xw_ar = join_rows(YLAG, Xw);
    if (which_loss==1){
      b = fastols(Y, Xw_ar, intercept);
    }
    if (which_loss==2){
      b = fastals(Y, Xw_ar, intercept, tau, 100, 1e-3);
    }
    if (which_loss==3){
      b = fastrq(Y, Xw_ar, intercept, tau);
    }
    if (intercept==1){
      Xw1 = join_rows(ones, Xw_ar);
    } else {
      Xw1 = Xw_ar;
    }
    res = Y-Xw1*b;
    if (which_loss==1){
      tmp = sum(pow(res,2))/n;
    }
    if (which_loss==2){
      arma::vec ind(n, fill::zeros);
      arma::vec ind2(n, fill::ones);
      arma::uvec idx = find(res < 0);
      ind(idx) = ind2(idx);
      tmp = sum(pow(res,2)%(tau - ind))/n;
    }
    if (which_loss==3){
      arma::vec ind(n, fill::zeros);
      arma::vec ind2(n, fill::ones);
      arma::uvec idx = find(res < 0);
      ind(idx) = ind2(idx);
      tmp = sum(res%(tau - ind))/n;
    }
    fit(i) = tmp(0,0);
  }
  arma::uvec idx2 = find(fit == min(fit));
  arma::vec par = param(idx2);
  w = pow(1-u,par(0)-1);
  w = w/sum(w);
  Xw = X*w;
  Xw_ar = join_rows(YLAG, Xw);
  if (which_loss==1){
    b = fastols(Y, Xw_ar, intercept);
  }
  if (which_loss==2){
    b = fastals(Y, Xw_ar, intercept, tau, 1000, 1e-7);
  }
  if (which_loss==3){
    b = fastrq(Y, Xw_ar, intercept, tau);
  }
  arma::vec beta = join_cols(b,par);
  
  return(beta);
}
