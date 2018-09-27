#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute the square-root Lasso
//'
//' Compute the square-root Lasso solution and its residuals for the nodewise Lasso.
//' @param y Response variable.
//' @param X Matrix of regressors.
//' @param lambda Regularization parameter.
//' @details The software is adapted from the Matlab-software provided by Belloni, Chernozhukov
//' and Wang (2011) in https://faculty.fuqua.duke.edu/~abn5/belloni-software.html
//' @export
// [[Rcpp::export]]
List sqrt_lasso(arma::vec y, arma::mat X, double lambda) {
  // Tolerance values
  double tolbeta = 1e-6;
  double tolobj = 1e-8;
  // Maximum number of iterations
  int maxiter = 10000;
  // Get dimensions of the regressor matrix
  int nObs = X.n_rows;
  int nVars = X.n_cols;
  // Initialize beta with the ridge regression estimator
  arma::mat RidgeMatrix = lambda * arma::eye<arma::mat>(nVars, nVars);
  arma::vec beta = arma::solve((X.t() * X) + RidgeMatrix, X.t() * y);
  arma::vec beta_old = beta;
  // Initialization
  int Iter = 0;
  // Gram matrix
  arma::mat XX = (X.t() * X) / nObs;
  arma::vec Xy = (X.t() * y) / nObs;
  arma::vec Error = y - (X * beta);
  double Qhat = pow(norm(Error, 2), 2) / nObs;
  double Snull = 0;
  double fobj = 0;
  double dual = 0;
  arma::vec oneVec(nVars, arma::fill::ones);
  arma::vec auxVec;
  while (Iter < maxiter) {
    Iter += 1;
    beta_old = beta;
    for (int j = 0; j < nVars; j++) {
      // Compute shoot
      Snull = arma::as_scalar((XX.row(j) * beta) - (XX(j, j) * beta(j)) - Xy(j));
      // Update beta_j
      if (fabs(beta(j)) > 0.0) {
        Error = Error + X.col(j) * beta(j);
        Qhat = pow(norm(Error, 2), 2) / nObs;
      }
      if (pow(nObs, 2) < pow(lambda, 2) / XX(j, j)) {
        beta(j) = 0.0;
      } else if (Snull > (lambda / nObs) * sqrt(Qhat)) {
        beta(j) = ( ( lambda / sqrt(pow(nObs, 2) - (pow(lambda, 2) / XX(j, j)))) *
          sqrt( std::max(Qhat - (pow(Snull, 2) / XX(j,j)), 0.0)) - Snull) / XX(j, j);
        Error = Error - X.col(j) * beta(j);
      } else if (Snull < - (lambda / nObs) * sqrt(Qhat)) {
        beta(j) = ( - ( lambda / sqrt(pow(nObs, 2) - (pow(lambda, 2) / XX(j, j)))) *
          sqrt( std::max(Qhat - (pow(Snull, 2) / XX(j,j)), 0.0)) - Snull) / XX(j, j);
        Error = Error - X.col(j) * beta(j);
      } else {
        beta(j) = 0.0;
      }
    }
    // Update primal and dual value
    // Primal
    fobj = norm((X * beta) - y, 2) / sqrt(nObs) + (lambda / nObs) * norm(beta, 1);
    // Dual
    if (norm(Error, 2) > 1e-10) {
      auxVec = (sqrt(nObs) * Error) / norm(Error, 2);
      dual = arma::as_scalar(((auxVec.t() * y) / nObs) - arma::abs(((lambda / nObs) * oneVec) -
        (arma::abs(X.t() * auxVec) / nObs)).t() * arma::abs(beta));
    } else {
      dual = (lambda / nObs) * norm(beta, 1);
    }
    // Stopping criterion
    if (norm(beta - beta_old, 1) < tolbeta) {
      if (fobj - dual < tolobj) {
        break;
      }
    }
  }
  arma::vec beta_out = beta;
  arma::vec Z = y - X * beta_out;
  return List::create(Named("beta") = beta_out,
                      Named("Z") = Z);
}
