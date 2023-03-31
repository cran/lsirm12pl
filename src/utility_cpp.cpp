// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
double Kumulavsech(arma::mat CC, double omega, arma::vec xlim, arma::vec ylim){
  double S = 0;
  int n = CC.n_rows;
  arma::vec s(n, fill::zeros);

  for(int i = 0; i < n; i++){
    s(i) = (normcdf(xlim(1), CC(i,0), omega) - normcdf(xlim(0), CC(i,0), omega))
    * (normcdf(ylim(1), CC(i,1), omega) - normcdf(ylim(0), CC(i,1), omega));
  }
  S = sum(s);

  return(S);
}

// [[Rcpp::export]]
double logpXCbeta(arma::mat X, arma::mat CC, double alpha, double omega, double AreaW, double integral){
  int nXrows = X.n_rows;
  int nCCrows = CC.n_rows;
  rowvec d;
  arma::vec res(nXrows, fill::zeros);
  arma::vec dummy;

  for(int i = 0; i < nXrows; i++){
    for(int j = 0; j < nCCrows; j++){
      d =  X.row(i) - CC.row(j);
      dummy =  (d*trans(d));
      res[i] = res[i] + exp(-dummy[0]/(2*omega*omega));
    }
  }
  double lik = AreaW - alpha*integral + nXrows*log(alpha/(2*(M_PI)*omega*omega)) + sum(log(res));

  return(lik);
}

// [[Rcpp::export]]
arma::vec NewPoint(arma::vec xlim, arma::vec ylim){
  double x, y = 0.0;
  arma::vec P(size(xlim), fill::zeros);

  x = randu(), y = randu();
  x = (xlim(1) - xlim(0)) * x + xlim(0);
  y = (ylim(1) - ylim(0)) * y + ylim(0);
  P(0) = x, P(1) = y;

  return(P);
}


