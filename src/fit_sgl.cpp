// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

typedef std::vector<int> stdint;
typedef std::vector<double> stdvec;

using namespace std;
//!--------------------------------- SGL solver -----------------------------------------!//
void linGradCalc(int *nrow, double *eta, double *y, double *ldot){
  for(int i = 0; i < nrow[0]; i++)
  {
    ldot[i] = (eta[i] - y[i])/nrow[0];  /* OR MAYBE NOT? */
  }
}

double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y){
  double squareSum = 0;
  
  for(int i = 0; i < nrow[0]; i++)
  {
    squareSum = squareSum + pow(eta[i] - y[i], 2)/2; 
  }
  
  return squareSum/nrow[0];   /* OR MAYBE NOT? */
}

// these functions allows for lambda1 and lambda2 to be vectors, i.e. penalty terms for each coefficients and group of cofficients

void linSolverLamGrid(double *X, double *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset){
  double *theta = new double[ncol[0]];
  int startInd = 0;
  double zeroCheck = 0;
  double check = 0;
  int count = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew = 0;
  double Lold = 0;
  double sqNormG = 0;
  double iProd = 0;
  double *etaNew = NULL;
  etaNew = new double[nrow[0]];
  double *etaNull = NULL;
  etaNull = new double[nrow[0]];
  //  int reset = 20;
  for(int i = 0; i < numGroup[0]; i++){ // loop through groups
    
    if(useGroup[i] == 1)
    {
      startInd = rangeGroupInd[i];
      
      // Setting up null gradient calc to check if group is 0
      for(int k = 0; k < nrow[0]; k++)
      {
        etaNull[k] = eta[k];
        for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
        {
          etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
        }
      }
      
      // Calculating Null Gradient
      linGradCalc(nrow, etaNull, y, ldot);
      
      double *grad = NULL;
      grad = new double[groupLen[i]];
      
      for(int j = 0; j < groupLen[i]; j++){ // within group loop
        grad[j] = 0;
        for(int k = 0; k < nrow[0]; k++)
        {
          grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
        }
        if(grad[j] < lambda1[j + rangeGroupInd[i]] && grad[j] > -lambda1[j + rangeGroupInd[i]])
        {
          grad[j] = 0;
        }
        if(grad[j] > lambda1[j + rangeGroupInd[i]])
        {
          grad[j] = grad[j] - lambda1[j + rangeGroupInd[i]];
        }
        if(grad[j] < -lambda1[j + rangeGroupInd[i]])
        {
          grad[j] = grad[j] + lambda1[j + rangeGroupInd[i]];
        }
        if(pow(grad[j],2) == pow(lambda1[j + rangeGroupInd[i]],2))
        {
          grad[j] = 0;
        }
      } // end of within group loop
      
      zeroCheck = 0;
      for(int j = 0; j < groupLen[i]; j++)
      {
        zeroCheck = zeroCheck + pow(grad[j],2);
      }
      
      if(zeroCheck <= pow(lambda2[i],2)*groupLen[i]){
        if(betaIsZero[i] == 0){
          for(int k = 0; k < nrow[0]; k++){
            for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++){
              eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
            }
          }
        }
        betaIsZero[i] = 1;
        for(int j = 0; j < groupLen[i]; j++){
          beta[j + rangeGroupInd[i]] = 0;
        }
      } else {
        if(isActive[i] == 0){
          groupChange = 1;
        }
        isActive[i] = 1;
        for(int k = 0; k < ncol[0]; k++){
          theta[k] = beta[k];
        }
        betaIsZero[i] = 0;
        double *z = NULL;
        z = new double[groupLen[i]];
        double *U = NULL;
        U = new double[groupLen[i]];
        double *G = NULL;
        G = new double[groupLen[i]];
        double *betaNew = NULL;
        betaNew = new double[ncol[0]];
        count = 0;
        check = 100000;
        while(count <= innerIter[0] && check > thresh[0]){
          count++;
          linGradCalc(nrow, eta, y ,ldot);
          for(int j = 0; j < groupLen[i]; j++){		  
            grad[j] = 0;
            for(int k = 0; k < nrow[0]; k++){
              grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
            }
          }
          diff = -1;
          Lold = linNegLogLikelihoodCalc(nrow, eta, y);
          
          // Back-tracking
          while(diff < 0){
            for(int j = 0; j < groupLen[i]; j++){
              z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
              if(z[j] < lambda1[j + rangeGroupInd[i]] * t && z[j] > -lambda1[j + rangeGroupInd[i]] * t){
                z[j] = 0;
              }
              if(z[j] > lambda1[j + rangeGroupInd[i]] * t){
                z[j] = z[j] - lambda1[j + rangeGroupInd[i]] * t;
              }
              if(z[j] < -lambda1[j + rangeGroupInd[i]] * t){
                z[j] = z[j] + lambda1[j + rangeGroupInd[i]] * t;
              }
            }
            norm = 0;
            for(int j = 0; j < groupLen[i]; j++){
              norm = norm + pow(z[j],2);
            }
            norm = sqrt(norm);
            if(norm != 0){
              uOp = (1 - lambda2[i]*sqrt(double(groupLen[i]))*t/norm);
            } else{
              uOp = 0;
            }
            if(uOp < 0){
              uOp = 0;
            }
            for(int j = 0; j < groupLen[i]; j++){
              U[j] = uOp*z[j];
              G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
            }
            
            // Setting up betaNew and etaNew in direction of Grad for descent step
            for(int k = 0; k < nrow[0]; k++){
              etaNew[k] = eta[k];
              for(int j = 0; j < groupLen[i]; j++){
                etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
              }
            }
            Lnew = linNegLogLikelihoodCalc(nrow, etaNew, y);
            sqNormG = 0;
            iProd = 0;
            for(int j = 0; j < groupLen[i]; j++){
              sqNormG = sqNormG + pow(G[j],2);
              iProd = iProd + grad[j] * G[j];
            }
            diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
            t = t * gamma[0];
          }
          t = t / gamma[0];
          check = 0;
          
          for(int j = 0; j < groupLen[i]; j++){
            check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
            for(int k = 0; k < nrow[0]; k++){
              eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
            }
            beta[j + rangeGroupInd[i]] = U[j] + count%reset[0]/(count%reset[0]+3) * (U[j] - theta[j + rangeGroupInd[i]]);
            theta[j + rangeGroupInd[i]] = U[j];
            
            for(int k = 0; k < nrow[0]; k++){
              eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
            }
          }
        }
        delete [] z;
        delete [] U;
        delete [] G;
        delete [] betaNew;
      }
      delete [] grad;
    }
  }
  delete [] etaNew;
  delete [] etaNull;
  delete [] theta;
}

void linNestLamGrid(double *X, double*y, int* index, int *nrow, int *ncol, 
                    int *numGroup, int *rangeGroupInd, int *groupLen, 
                    double *lambda1, double *lambda2, double *beta, 
                    int *innerIter, int *outerIter, double *thresh, 
                    double *outerThresh, double *eta, double *gamma, 
                    int *betaIsZero, double *step, int *reset){
  // depends on: linSolver;
  double* prob = NULL;
  prob = new double[nrow[0]];
  double* nullBeta = NULL;
  nullBeta = new double[ncol[0]];
  int n = nrow[0];
  int p = ncol[0];
  double *ldot = NULL;
  ldot = new double[n];
  int groupChange = 1;
  int* isActive = NULL;
  isActive = new int[numGroup[0]];
  int* useGroup = NULL;
  useGroup = new int[numGroup[0]];
  int* tempIsActive = NULL;
  tempIsActive = new int[numGroup[0]];
  
  for(int i = 0; i < numGroup[0]; i++){
    isActive[i] = 0;
    useGroup[i] = 1;
  }
  
  // outer most loop creating response etc...
  int outermostCounter = 0;
  double outermostCheck = 100000;
  double* outerOldBeta = NULL;
  outerOldBeta = new double[p];
  
  while(groupChange == 1){
    groupChange = 0;
    linSolverLamGrid(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, step, reset);
    while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0]){
      outermostCounter ++;
      for(int i = 0; i < p; i++){
        outerOldBeta[i] = beta[i];
      }
      for(int i = 0; i < numGroup[0]; i++){
        tempIsActive[i] = isActive[i];
      }
      
      linSolverLamGrid(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, step, reset);
      
      outermostCheck = 0;
      for(int i = 0; i < p; i++){
        outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
      }
    }
  }
  
  delete [] nullBeta;
  delete [] outerOldBeta;
  delete [] ldot;
  delete [] isActive;
  delete [] useGroup;
  delete [] tempIsActive;
  delete [] prob;
  
  //return 1;
}

//!--------------------------------- Other functions -----------------------------------------!//
//' Computes cross-validation folds. 
//' 
//' @param nfolds (double) number of folds.
//' @param nrow (double) number of rows of data. 
//' @return foldsid (vec) vector of cv folds.
//' @export cvfolds
//' @keywords internal
// [[Rcpp::export]]
arma::vec cvfolds(double nfolds, double nrow){
  mat tmp_foldsid = linspace(0,nfolds-1,nfolds);
  mat tmp2_foldsid = repmat(tmp_foldsid, 1, ceil(nrow/nfolds));
  vec foldsid = vectorise(tmp2_foldsid);
  foldsid = foldsid.subvec(0,nrow-1);
  foldsid = sort(foldsid);
  return foldsid;
}

//' Computes optimal \eqn{\lambda} for cross-validation output.
//' 
//' @param lambda (vec) vector of \eqn{\lambda} parameters.
//' @param cvm (vec) vector of cv error curve. 
//' @param cvsd (vec) vector of standard errors for cv error curve.
//' @param which_lambda which optimal \eqn{\lambda} to select; 0 - min of cv error curve, 1 - 1se min of cv error curve.
//' @export getmin_cpp
//' @keywords internal
// [[Rcpp::export]]
double getmin_cpp(arma::vec lambda,arma::vec cvm,arma::vec cvsd, int which_lambda){
  // which_lambda = 0 - CV-curve minimizing lambda; 
  // which_lambda = 1 - CV-curve+1std minimizing lambda;
  double cvmin=min(cvm);
  double lambda_min = lambda(max(find(cvm<=cvmin)));
  arma::vec cvm_plus_cvsd = cvm+cvsd;
  double lambda_1se = lambda(max(find(cvm<=cvm_plus_cvsd(max(find(cvm<=cvmin))))));
  double lambda_opt=0;
  if (which_lambda==0){
    lambda_opt = lambda_min;
  } 
  if (which_lambda==1){
    lambda_opt = lambda_1se;
  }
  if (which_lambda>2){
    lambda_opt = datum::inf;
  }
  return lambda_opt;
}

//' Computes a single solution for \eqn{\lambda_1} and \eqn{\lambda_2} values
//' 
//' @param beta0 (vec) vector of initial \eqn{\beta_0} regression coefficients.
//' @param Z (mat) matrix of dummies. 
//' @param X (mat) matrix of covariates.
//' @param y (vec) vector of response.
//' @param index (vec) vector indicating group membership of each covariate.
//' @param lambda1 (vec) single value of \eqn{\lambda_1}.
//' @param lambda2 (vec) single value of \eqn{\lambda_2}.
//' @param innerIter (int) max number of inner iterations.
//' @param outerIter (int) max number of outer iterations.
//' @param thresh (double) convergence threshold of inner loop.
//' @param outerThresh (double) convergence threshold of outer loop.
//' @param gamma_solver (double) solver parameter.
//' @param step (double) solver parameter.
//' @param reset (int) solver parameter.
//' @export cpp_sgl_fit
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_sgl_fit(arma::vec& beta0, arma::mat& Z, arma::mat& X, arma::vec& y, arma::vec& index,
                      arma::vec& lambda1, arma::vec& lambda2, 
                      int innerIter, int outerIter, double thresh, 
                      double outerThresh, double gamma_solver, 
                      double step, int reset)  { 
  // ====================================================================
  // purpose: wrap up linNest C++ function for minimizing mse-sg-LASSO:
  // 
  //   min_(alpha,beta) ||Zalpha + Xbeta - y||^2_T + lambda ( gamma_w*|beta|_1 + (1-gamma_w)*|beta|_1,2 )
  //
  // using gradient/block-wise coordinate descent of SGL R package.
  // This function computes the full path of lambdas, 
  // for a fixed gamma_w weight (lasso vs group relative weight)
  // lambda1 and lambda2 are vectors of size p and |G| respectively, 
  // where |G| denotes the number of groups.
  // ZX - matrix of dummies and covariates
  // ====================================================================
  double ndum = Z.n_cols;
  arma::mat  X1 = join_rows(Z, X);
  int ncol = X1.n_cols, nrow = X1.n_rows;
  
  // sg-LASSO stuff:
  arma::vec groups = unique(index);
  int numGroup = groups.n_elem;
  arma::vec rangeGroupInd(numGroup+1, fill::zeros);
  for (int i = 0; i < numGroup; i ++){
    rangeGroupInd[i] =  min(find(index == groups[i]));
  }
  rangeGroupInd[numGroup] = ncol;
  arma::vec groupLen = diff(rangeGroupInd);
  
  arma::vec betaIsZero(numGroup,fill::ones);
  arma::vec eta(nrow,fill::zeros);
  arma::vec beta_int = beta0;//(ncol,fill::zeros);;
  
  // convert variables to format that fits solver:
  arma::vec Xvec = vectorise(X1);
  stdvec pass_X = arma::conv_to< stdvec >::from(Xvec);  
  stdvec pass_y = arma::conv_to< stdvec >::from(y);
  stdint pass_index = arma::conv_to< stdint >::from(index);
  stdint pass_rangeGroupInd = arma::conv_to< stdint >::from(rangeGroupInd);
  stdint pass_groupLen = arma::conv_to< stdint >::from(groupLen);
  stdint pass_betaIsZero = arma::conv_to< stdint >::from(betaIsZero);
  stdvec pass_eta = arma::conv_to< stdvec >::from(eta);
  stdvec pass_beta = arma::conv_to< stdvec >::from(beta_int);
  
  double *p_X = &pass_X[0];
  double *p_y = &pass_y[0];
  int *p_index = &pass_index[0];
  int *p_rangeGroupInd = &pass_rangeGroupInd[0];
  int *p_groupLen = &pass_groupLen[0];
  double *p_eta = &pass_eta[0];
  int *p_betaIsZero = &pass_betaIsZero[0];
  double *p_beta = &pass_beta[0];

  // solve sg-LASSO path
  stdvec pass_lambda1 = arma::conv_to< stdvec >::from(lambda1);
  stdvec pass_lambda2 = arma::conv_to< stdvec >::from(lambda2);
  double *p_lambda1 = &pass_lambda1[0];
  double *p_lambda2 = &pass_lambda2[0];
  // solve MSE sg-LASSO using gradient descent
  linNestLamGrid(p_X, p_y, p_index, &nrow, &ncol, &numGroup, p_rangeGroupInd, p_groupLen, p_lambda1, p_lambda2, p_beta, &innerIter, &outerIter, &thresh, &outerThresh, p_eta, &gamma_solver, p_betaIsZero, &step, &reset);
  // store solution
  arma::vec out_beta(ncol);
  for (int i = 0; i < ncol; i ++) {
    out_beta(i) = p_beta[i];    
  }
  if (ndum==1){
      double tmp_m = mean(y);
      out_beta.subvec(0,ndum-1) = out_beta.subvec(0,ndum-1) - tmp_m;
  } else {
    for (int i = 0; i < ndum; i++){
      double tmp_m = mean(y%Z.col(i));
      out_beta.subvec(i,i) = out_beta.subvec(i,i) - tmp_m;
    }
  }
 return out_beta;
}

//' Computes a path solution for \eqn{\lambda} and fixed \eqn{\gamma} values
//' 
//' @param X (mat) matrix of covariates.
//' @param Z (mat) matrix of dummies. 
//' @param y (vec) vector of response.
//' @param index (vec) vector indicating group membership of each covariate.
//' @param dummies (double) add dummy variables in a regression (dummies = 1) or not (dummies = 0).
//' @param l1_frac \eqn{\ell_1} norm peanlty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_1 norm}).
//' @param l21_frac \eqn{\ell_{1,2}} norm peanlty factor for random or fixed effects (default value is zero which means \eqn{\alpha} is left unpenalized in \eqn{\ell_{1,2}} norm).
//' @param dummies_index vector indicating group membership of \eqn{\alpha}.
//' @param lambdas (vec) sequence of \eqn{\lambda} values for fitting.
//' @param gamma_w (double) sg-LASSO mixing parameter.
//' @param innerIter (int) max number of inner iterations.
//' @param outerIter (int) max number of outer iterations.
//' @param thresh (double) convergence threshold of inner loop.
//' @param outerThresh (double) convergence threshold of outer loop.
//' @param gamma_solver (double) solver parameter.
//' @param step (double) solver parameter.
//' @param reset (int) solver parameter.
//' @export cpp_sgl_fitpath
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_sgl_fitpath(arma::mat& X, arma::mat& Z, arma::vec& y, arma::vec& index, double dummies,
                      arma::vec l1_frac, arma::vec l21_frac, arma::vec dummies_index,
                      arma::vec& lambdas, double gamma_w, 
                      int innerIter, int outerIter, double thresh, 
                      double outerThresh,double gamma_solver, 
                      double step, int reset)  { 
  // arma::vec zero_pen(ndum, fill::zeros);
  // len_lam1 = join_cols(zero_pen,len_lam1);
  // len_lam2 = join_cols(zero_pen,len_lam2);
  // arma::vec dummies_index = linspace(1,ndum,ndum);
  // index = join_cols(dummies_index,index + ndum);
  // ====================================================================
  // purpose: wrap up linNest C++ function for minimizing mse-sg-LASSO:
  // 
  //   min_x ||Ax - b||^2_T + lambda1 * |x|_1 + lambda2* |x|_1,2
  // using gradient/block-wise coordinate descent of SGL R package.
  // This function computes the full path of lambdas, 
  // for a fixed gamma_w weight (lasso vs group relative weight)
  // here lambda1 and lambda2 are vectors of size p and |G| respectively, 
  // where |G| denotes the number of groups.
  // ====================================================================
  int nrow = X.n_rows;
  arma::vec len_lam1(X.n_cols,fill::ones);
  arma::vec len_lam2(max(index),fill::ones);
  arma::mat X1;
  double ndum = Z.n_cols;
  if (dummies == 1){
    len_lam1 = join_cols(l1_frac,len_lam1);
    len_lam2 = join_cols(l21_frac,len_lam2);
    index = join_cols(dummies_index,index + max(dummies_index));
    X1 = join_rows(Z, X);
  } else {
    X1 = X;
  }
  int ncol = X1.n_cols, nlam = lambdas.n_elem;
  
  // sg-LASSO stuff:
  arma::vec groups = unique(index);
  int numGroup = groups.n_elem;
  arma::vec rangeGroupInd(numGroup+1, fill::zeros);
  for (int i = 0; i < numGroup; i ++){
    rangeGroupInd[i] =  min(find(index == groups[i]));
  }
  rangeGroupInd[numGroup] = ncol;
  arma::vec groupLen = diff(rangeGroupInd);
  
  arma::vec betaIsZero(numGroup,fill::ones);
  arma::vec eta(nrow,fill::zeros);
  arma::vec beta_int(ncol, fill::zeros);
  if(dummies==1){
    double nd = Z.n_cols;
    beta_int.subvec(0,nd-1)=Z.t()*y/nrow;
  }
  
  // convert variables to format that fits solver:
  arma::vec Xvec = vectorise(X1);
  stdvec pass_X = arma::conv_to< stdvec >::from(Xvec);  
  stdvec pass_y = arma::conv_to< stdvec >::from(y);
  stdint pass_index = arma::conv_to< stdint >::from(index);
  stdint pass_rangeGroupInd = arma::conv_to< stdint >::from(rangeGroupInd);
  stdint pass_groupLen = arma::conv_to< stdint >::from(groupLen);
  stdint pass_betaIsZero = arma::conv_to< stdint >::from(betaIsZero);
  stdvec pass_eta = arma::conv_to< stdvec >::from(eta);
  stdvec pass_beta = arma::conv_to< stdvec >::from(beta_int);
  
  double *p_X = &pass_X[0];
  double *p_y = &pass_y[0];
  int *p_index = &pass_index[0];
  int *p_rangeGroupInd = &pass_rangeGroupInd[0];
  int *p_groupLen = &pass_groupLen[0];
  double *p_eta = &pass_eta[0];
  int *p_betaIsZero = &pass_betaIsZero[0];
  double *p_beta = &pass_beta[0];
  
  arma::mat out_betas(ncol,nlam);
  // solve sg-LASSO path
  for (int k = 0; k < nlam; k ++){
    Rcpp::checkUserInterrupt();
    arma::vec lambda1 = len_lam1*lambdas(k)*gamma_w;
    arma::vec lambda2 = len_lam2*lambdas(k)*(1-gamma_w);
    stdvec pass_lambda1 = arma::conv_to< stdvec >::from(lambda1);
    stdvec pass_lambda2 = arma::conv_to< stdvec >::from(lambda2);
    double *p_lambda1 = &pass_lambda1[0];
    double *p_lambda2 = &pass_lambda2[0];
    // solve MSE sg-LASSO using gradient descent
    linNestLamGrid(p_X, p_y, p_index, &nrow, &ncol, &numGroup, p_rangeGroupInd, p_groupLen, p_lambda1, p_lambda2, p_beta, &innerIter, &outerIter, &thresh, &outerThresh, p_eta, &gamma_solver, p_betaIsZero, &step, &reset);
    // store kth solution
    arma::vec out_beta(ncol);
    for (int i = 0; i < ncol; i ++) {
      out_beta(i) = p_beta[i];
    }
    if (dummies==1){
      if (ndum==1){
        double tmp_m = mean(y);
        out_beta.subvec(0,ndum-1) = out_beta.subvec(0,ndum-1) - tmp_m;
      } else {
        for (int i = 0; i < ndum; i++){
          double tmp_m = mean(y%Z.col(i));
          out_beta.subvec(i,i) = out_beta.subvec(i,i) - tmp_m;
        }
      }
    }
    out_betas.col(k) = out_beta;
  }
  return out_betas;
}