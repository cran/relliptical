// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>

#include <Rdefines.h>

using namespace Rcpp;
using namespace arma;
using namespace Numer;


// RANDOM NUMBER GENERATOR

// Generate random numbers from Truncated Multivariate Student's-t distribution t(mu,Sigma,nu; (a,b))
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat rttrunc(arma::uword n, double nu, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  double nup = nu + p;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s * s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  Rcpp::NumericVector l1 = Rcpp::wrap((a - mu)/s);
  Rcpp::NumericVector u1 = Rcpp::wrap((b - mu)/s);
  arma::vec pa = Rcpp::pt(l1,nu,1,0);
  arma::vec pb = Rcpp::pt(u1,nu,1,0);
  arma::vec x0 = randu<arma::vec>(p);
  Rcpp::NumericVector x1 = Rcpp::wrap(pa + (pb - pa)%x0);
  arma::colvec x = Rcpp::qt(x1,nu,1,0);
  arma::vec lower = Rcpp::as<arma::vec>(l1);
  arma::vec upper = Rcpp::as<arma::vec>(u1);

  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = lower.elem(q1);
  q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    y = arma::randu<double>()*std::exp(-0.5*nup*std::log(1.0+delta/nu));
    kap = nu*(std::pow(y,-2.0/nup) - 1.0);
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from Truncated Multivariate Normal distribution N(mu,Sigma; (a,b))
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat rtnormal(arma::uword n, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s*s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  Rcpp::NumericVector l1 = Rcpp::wrap((a - mu)/s);
  Rcpp::NumericVector u1 = Rcpp::wrap((b - mu)/s);
  arma::vec pa = Rcpp::pnorm(l1,0,1,1,0);
  arma::vec pb = Rcpp::pnorm(u1,0,1,1,0);
  arma::vec x0 = randu<arma::vec>(p);
  Rcpp::NumericVector x1 = Rcpp::wrap(pa + (pb - pa)%x0);
  arma::colvec x = Rcpp::qnorm(x1,0,1,1,0);
  arma::vec lower = Rcpp::as<arma::vec>(l1);
  arma::vec upper = Rcpp::as<arma::vec>(u1);

  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = lower.elem(q1);
  q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    kap = -2.0*std::log(arma::randu<double>()) + delta;
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj + (kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv - lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from Truncated Multivariate Power Exponential distribution PE(mu,Sigma,beta; (a,b))
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat rtPE(arma::uword n, double beta, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s*s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  arma::vec lower = (a - mu)/s;
  arma::vec upper = (b - mu)/s;
  arma::vec x = lower;
  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);
  x.replace(arma::datum::inf,0.0);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    y = -2.0*std::log(arma::randu<double>()) + std::pow(delta,beta);
    kap = std::pow(y, 1.0/beta);
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from Truncated Multivariate Pearson VII distribution PVII(mu,Sigma,N,nu; (a,b))
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat rtPVII(arma::uword n, double N, double nu, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s*s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  Rcpp::NumericVector l1 = Rcpp::wrap((a - mu)/s);
  Rcpp::NumericVector u1 = Rcpp::wrap((b - mu)/s);
  double Ni = N - (p-1)/2.0;
  double nui = sqrt((2.0*Ni-1.0)/nu);
  arma::vec pa = Rcpp::pt(nui*l1,2.0*Ni-1.0,1,0);
  arma::vec pb = Rcpp::pt(nui*u1,2.0*Ni-1.0,1,0);
  arma::vec x0 = randu<arma::vec>(p);
  Rcpp::NumericVector x1 = Rcpp::wrap(pa + (pb - pa)%x0);
  arma::colvec x = (Rcpp::qt(x1,2.0*Ni-1.0,1,0))/nui;
  arma::vec lower = Rcpp::as<arma::vec>(l1);
  arma::vec upper = Rcpp::as<arma::vec>(u1);

  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = lower.elem(q1);
  q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    y = arma::randu<double>()*std::exp(-N*std::log(1.0 + delta/nu));
    kap = nu*(std::pow(y,-1.0/N) - 1.0);
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from Truncated Multivariate Slash distribution Slash(mu,Sigma,nu; (a,b))
// ------------------------------------------------------------------------------------------------------------
class Mintegrand: public Func
{
private:
  const double q; // degrees of freedom nu
  const double p; // length of multivariate vector p
  const double u; // Mahalanobis distance
public:
  Mintegrand(double q_, double p_, double u_) : q(q_), p(p_), u(u_) {}
  double operator()(const double& x) const
  {
    return std::pow(x, (q+0.50*p-1.0))*std::exp(-0.50*u*x);
  }
};

// Compute the slash DGF for a given values (u,nu,p)
double slash_g(double u, double nu, double p){
  Mintegrand f(nu, p, u);
  double err_est;
  int err_code;
  return integrate(f, 0.0, 1.0, err_est, err_code);
}

// Compute the inverse of the slash DGF using Brent's algorithm
double r8_epsilon ( ){
  const double value = 2.220446049250313E-016;
  return value;
}
double BrentMethod(double y1, double nu, double p1){
  double a = 0.0001; double b = 1e5; double t = 1e-10;
  double c; double d; double e; double fa; double fb; double fc;
  double m; double tol; double macheps = r8_epsilon( );
  double p; double q; double r; double s; double sa; double sb;

  sa = a;
  sb = b;
  fa = slash_g(sa,nu,p1) - y1;
  fb = slash_g(sb,nu,p1) - y1;
  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;
  if (fa * fb > 0.0) stop("f() values at end points not of opposite sign.");

  for ( ; ; )
  {
    if (abs(fc) < abs(fb)){
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol = 2.0*macheps*abs(sb) + t;
    m = 0.5*(c - sb);

    if (abs(m) <= tol || fb == 0.0 ){
      break;
    }

    if (abs(e) < tol || abs(fa) <= abs(fb)){
      e = m; d = e;
    } else {
      s = fb/fa;
      if (sa == c){
        p = 2.0*m*s;
        q = 1.0 - s;
      } else {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*m*q*(q - r) - (sb - sa)*(r - 1.0));
        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      }
      if ( 0.0 < p ){
        q = - q;
      } else {
        p = - p;
      }
      s = e;
      e = d;
      if (2.0*p < 3.0*m*q - abs(tol*q) && p < abs(0.5*s*q)){
        d = p/q;
      } else {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;
    if (tol < abs(d)){
      sb = sb + d;
    } else if ( 0.0 < m ){
      sb = sb + tol;
    } else {
      sb = sb - tol;
    }
    fb = slash_g(sb,nu,p1) - y1;
    if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)){
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}

// [[Rcpp::export]]
arma::mat rtslash(arma::uword n, double nu, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s*s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  arma::vec lower = (a - mu)/s;
  arma::vec upper = (b - mu)/s;
  arma::vec x = lower;
  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);
  x.replace(arma::datum::inf,0.0);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    y = arma::randu<double>()*slash_g(delta,nu,p);
    kap = BrentMethod(y,nu,p);
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from Truncated Multivariate Contaminated Normal distribution CN(mu,Sigma,nu,rho; (a,b))
// ----------------------------------------------------------------------------------------------------------------

// Compute the inverse for the Contaminated normal DGF using Newton Raphson
double ginvCN(double nu, double rho, double p, double y){
  double t1 = -2.0*std::log(y);
  double c1 = nu*std::pow(rho,(0.5*p));
  double t2 = t1 + 2.0*(c1*std::exp(-0.5*rho*t1) + (1.0-nu)*std::exp(-0.5*t1) - y)/(rho*c1*std::exp(-0.5*rho*t1) + (1.0-nu)*std::exp(-0.5*t1));
  while (abs(t2 - t1) > 1e-10){
    t1 = t2;
    t2 = t1 + 2.0*(c1*std::exp(-0.5*rho*t1) + (1.0-nu)*std::exp(-0.5*t1) - y)/(rho*c1*std::exp(-0.5*rho*t1) + (1.0-nu)*std::exp(-0.5*t1));
  }
  if (t2 < 0.0){ t2 = 0.0;}
  return t2;
}

// [[Rcpp::export]]
arma::mat rtCN(arma::uword n, double nu, double rho, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, int burn, int lag){
  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s*s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  arma::vec lower = (a - mu)/s;
  arma::vec upper = (b - mu)/s;
  arma::vec x = lower;
  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);
  x.replace(arma::datum::inf,0.0);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, dgf, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    dgf = nu*std::pow(rho,(0.5*p))*std::exp(-0.5*rho*delta) + (1.0-nu)*std::exp(-0.5*delta);
    y = arma::randu<double>()*(dgf);
    kap = ginvCN(nu,rho,p,y);
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Generate random numbers from any Truncated Multivariate distribution given DGF f(mu,Sigma; gFun, (a,b))
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat randomG(arma::uword n, arma::vec mu, arma::mat Sigma, arma::vec a, arma::vec b, Function gFUN, Function ginvFUN, int burn, int lag){

  arma::uword m = lag*n + burn;
  arma::uword p = Sigma.n_cols;
  arma::vec s = sqrt(Sigma.diag());
  arma::mat R = Sigma%(1.0/(s * s.t()));
  arma::mat Rinv = R.i();
  arma::mat X(n,p,fill::zeros);

  arma::vec lower = (a - mu)/s;
  arma::vec upper = (b - mu)/s;
  arma::vec x = lower;
  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);
  x.replace(arma::datum::inf,0.0);

  arma::umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, y, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  arma::uword count = 1;
  for(uword i=0;i<m;i++){
    delta = as_scalar(x.t()*Rinv*x);
    y = arma::randu<double>()*as<double>(gFUN(delta));
    kap = as<double>(ginvFUN(y));
    for(uword j=0;j<p;j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj+(kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv-lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// COMPUTE MOMENTS - E(Y), E(YYt), Var(Y)

// Compute moments from truncated multivariate Student-t distribution
// ------------------------------------------------------------------------------------------------------------
class Myintegral: public Func
{
private:
  const double c; // nu
public:
  Myintegral(double c_) : c(c_) {}
  double operator()(const double& x) const
  {
    return std::pow(1.0 + x*x, c);
  }
};

double E2integ(double nu, double a, double b){
  Myintegral g(nu);
  double err_est;
  int err_code;
  return integrate(g, a, b, err_est, err_code);
}

List tuniv(double mu, double sigma, double nu, arma::vec lower, arma::vec upper, arma::uword n2){
  double md0 = arma::datum::nan;
  double vr0 = arma::datum::nan;

  if (n2 == 1){ // Exists mean and variance
    double s11 = sqrt(sigma);
    double a = (as_scalar(lower) - mu)/s11;
    double b = (as_scalar(upper) - mu)/s11;
    double alpha0 = R::pt(b,nu,1,0) - R::pt(a,nu,1,0);
    double meanx = 0.0, varx = 0.0;
    if (nu == 1.0){
      meanx = log((1.0 + b*b)/(1.0 + a*a))/(2.0*M_PI*alpha0);
    } else {
      meanx = tgamma(0.5*(nu + 1.0))/(alpha0*sqrt(nu*M_PI)*tgamma(0.5*nu))*(nu/(1.0-nu))*(pow(1.0 + b*b/nu, -0.5*(nu-1.0)) - pow(1.0 + a*a/nu, -0.5*(nu-1.0)));
    }
    md0 = mu + s11*meanx; // Mean of Y
    if (nu > 2.0){
      varx = (nu*(nu-1.0)/((nu-2.0)*alpha0))*(R::pt(b*sqrt((nu-2.0)/nu),nu-2.0,1,0) - R::pt(a*sqrt((nu-2.0)/nu),nu-2.0,1,0)) - nu - meanx*meanx;
    } else {
      if (nu == 2.0){
        varx = (asinh(b/sqrt(2.0)) - asinh(a/sqrt(2.0)))/alpha0 - 2.0 - meanx*meanx;
      } else {
        varx = nu*sqrt(nu)*tgamma(0.5*(nu + 1.0))/(alpha0*sqrt(nu*M_PI)*tgamma(0.5*nu))*E2integ(-0.50*(nu-1.0),a/sqrt(nu),b/sqrt(nu)) - nu - meanx*meanx;
      }
    }
    vr0 = varx*(sigma);  // Variance  of Y
  } else {

    if (nu > 1.0){ // Exists the mean
      double s11 = sqrt(sigma);
      double a, b;
      if(lower.is_finite()){ a = (as_scalar(lower) - mu)/s11; } else { a = -1e40; }
      if(upper.is_finite()){ b = (as_scalar(upper) - mu)/s11; } else { b = 1e40; }
      double alpha0 = R::pt(b,nu,1,0) - R::pt(a,nu,1,0);
      double meanx = tgamma(0.5*(nu + 1.0))/(alpha0*sqrt(nu*M_PI)*tgamma(0.5*nu))*(nu/(1.0-nu))*(pow(1.0 + b*b/nu, -0.5*(nu-1.0)) - pow(1.0 + a*a/nu, -0.5*(nu-1.0)));
      md0 = mu + s11*meanx;
      if (nu > 2.0){ // Exists the variance
        double varx = (nu*(nu-1.0)/((nu-2.0)*alpha0))*(R::pt(b*sqrt((nu-2.0)/nu),nu-2.0,1,0) - R::pt(a*sqrt((nu-2.0)/nu),nu-2.0,1,0)) - nu - meanx*meanx;
        vr0 = varx*(sigma);
      }
    }
  }
  List out;
  out["mean"] = md0; out["var"] = vr0;
  return out;
}

// [[Rcpp::export]]
List Tmoment(arma::vec mu, arma::mat Sigma, double nu, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::vec mean0(p,fill::zeros);   mean0.replace(0,arma::datum::nan);
  arma::mat var0(p,p,fill::zeros);  var0.replace(0,arma::datum::nan);
  arma::mat mom20(p,p,fill::zeros); mom20.replace(0,arma::datum::nan);
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();

  if (lind == p){ // Non-truncated variables
    if (nu > 1.0){
      mean0 = mu;
      if (nu > 2.0){ var0 = Sigma*(nu)/(nu - 2.0); mom20 = var0 + mean0*mean0.t(); }
    }
  } else {

    if ((lind==0) & (p==1)){ // All variables are truncated: univariate case
      arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
      List momentos = tuniv(as_scalar(mu),as_scalar(Sigma),nu,lower,upper,ind0.size());
      mean0(0) = as<double>(momentos["mean"]);
      var0(0,0) = as<double>(momentos["var"]);
      mom20(0,0) = var0(0,0) + mean0(0)*mean0(0);
    } else {

      if ((lind==0) & (p>1)){ // All variables are truncated: p-variate case
        arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
        arma::uword n2 = ind0.size();
        double d = n2 + nu;

        if (d > 2.0){ // Exists mean and variance
          arma::mat gen = rttrunc(n,nu,mu,Sigma,lower,upper,burn,thinning);
          mean0 = (mean(gen,0)).t();
          var0 = cov(gen); mom20 = var0 + mean0*mean0.t();
        } else {
          if ((d>1.0) & (d<=2.0)){
            if (n2 == 0){ // Exists mean
              arma::mat gen = rttrunc(n,nu,mu,Sigma,lower,upper,burn,thinning);
              mean0 = (mean(gen,0)).t();
            } else { // Exists mean (all), variance (bounded), covariance (all)
              arma::uvec ind3 = unique(join_vert(find_nonfinite(lower),find_nonfinite(upper))); // Non-bounded region
              arma::mat Anan(ind3.size(),ind3.size(),fill::zeros); Anan.replace(0,arma::datum::nan);
              arma::mat gen = rttrunc(n,nu,mu,Sigma,lower,upper,burn,thinning);
              mean0 = (mean(gen,0)).t();
              var0 = cov(gen);                    var0(ind3,ind3) = Anan;
              mom20 = var0 + mean0*mean0.t();     mom20(ind3,ind3) = Anan;
            }
          }
        }
      } else {
        if ((lind==(p-1)) & (p>1)){ // One variable is truncated: p-variate case
          arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
          arma::uword n2 = ind0.size();
          double d = n2 + nu;

          if (d > 2.0){ // Exists mean and variance
            arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
            List momentos = tuniv(as_scalar(mu(ind1)),as_scalar(Sigma(ind1,ind1)),nu,lower(ind1),upper(ind1),n2);
            arma::vec mx(1);    mx(0) = as<double>(momentos["mean"]);
            mean0(ind1) = mx;
            mean0(ind2) = mu(ind2) + as_scalar((mean0(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
            arma::mat vx(1,1); vx(0,0) = as<double>(momentos["var"]);
            var0(ind1,ind1) = vx;
            double omega21 = (nu + as_scalar((var0(ind1,ind1) +  pow(mean0(ind1)-mu(ind1),2.0))/Sigma(ind1,ind1)))/(nu-1.0);
            var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - (omega21 - as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1)))/as_scalar(Sigma(ind1,ind1))*Sigma(ind2,ind1)*Sigma(ind1,ind2);
            var0(ind2,ind1) = as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1))*Sigma(ind2,ind1);
            var0(ind1,ind2) = (var0(ind2,ind1)).t();
            mom20 = var0 + mean0*mean0.t();
          } else {
            if ((d>1.0) & (d<=2.0)){
              if (n2 == 0){ // Exists mean
                arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
                List momentos = tuniv(as_scalar(mu(ind1)),as_scalar(Sigma(ind1,ind1)),nu,lower(ind1),upper(ind1),n2);
                arma::vec meanx(1); meanx(0) = as<double>(momentos["mean"]);
                mean0(ind1) = meanx;
                mean0(ind2) = mu(ind2) + as_scalar((mean0(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
              } else{ // Exists mean (all), variance (bounded), covariance (all)
                List momentos = tuniv(as_scalar(mu(ind0)),as_scalar(Sigma(ind0,ind0)),nu,lower(ind0),upper(ind0),n2);
                arma::vec meanx(1);  meanx(0) = as<double>(momentos["mean"]);
                mean0(ind0) = meanx;
                mean0(ind2) = mu(ind2) + as_scalar((mean0(ind0) - mu(ind0))/Sigma(ind0,ind0))*Sigma(ind2,ind0);
                arma::mat vx(1,1); vx(0,0) = as<double>(momentos["var"]);
                var0(ind0,ind0) = vx;
                var0(ind2,ind0) = as_scalar(var0(ind0,ind0)/Sigma(ind0,ind0))*Sigma(ind2,ind0);
                var0(ind0,ind2) = (var0(ind2,ind0)).t();
                mom20 = var0 + mean0*mean0.t();
              }
            }
          }
        } else { // The number of truncated variables varies between 2 and p-1
          arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper));         // Doubly truncated variables
          arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
          arma::uword n2 = ind0.size();
          double d = n2 + nu;
          if (d > 2.0){ // Exists mean and variance
            arma::mat sigInv = (Sigma(ind1,ind1)).i();
            arma::mat gen = rttrunc(n,nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
            mean0(ind1) = (mean(gen,0)).t();
            mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
            arma::mat Iden = eye(ind1.size(),ind1.size());
            var0(ind1,ind1) = cov(gen);
            double omega21 = (nu + trace(var0(ind1,ind1)*sigInv) + as_scalar((mean0(ind1) - mu(ind1)).t()*sigInv*(mean0(ind1) - mu(ind1))))/(nu + ind1.size() - 2.0);
            var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - Sigma(ind2,ind1)*sigInv*(omega21*Iden - var0(ind1,ind1)*sigInv)*Sigma(ind1,ind2);
            var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
            var0(ind1,ind2) = (var0(ind2,ind1)).t();
            mom20 = var0 + mean0*mean0.t();
          } else {
            if ((d>1.0) & (d<=2.0)){
              if (n2 == 0){ // Exists mean
                arma::mat sigInv = (Sigma(ind1,ind1)).i();
                arma::mat gen = rttrunc(n,nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
                mean0(ind1) = (mean(gen,0)).t();
                mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
              } else { // Exists mean (all), variance (bounded), covariance (all)
                arma::uvec ind3 = unique(join_vert(find_nonfinite(lower),find_nonfinite(upper))); // Non-bounded region
                arma::mat Anan(ind3.size(),ind3.size(),fill::zeros); Anan.replace(0,arma::datum::nan);
                arma::mat sigInv = (Sigma(ind1,ind1)).i();
                arma::mat gen = rttrunc(n,nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
                mean0(ind1) = (mean(gen,0)).t();
                mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
                var0(ind1,ind1) = cov(gen);
                var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
                var0(ind1,ind2) = (var0(ind2,ind1)).t();
                var0(ind3,ind3) = Anan;
                mom20 = var0 + mean0*mean0.t(); mom20(ind3,ind3) = Anan;
              }
            }
          }
        } // End if
      } // End if
    } // End if
  } // End if
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.5*(mom20 + mom20.t());
  output["VarY"] = 0.5*(var0 + var0.t());
  return output;
}


// Compute moments from truncated multivariate Normal distribution
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List Nmoment(arma::vec mu, arma::mat Sigma, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();
  arma::vec mean0(p); arma::mat var0(p,p); arma::mat mom20(p,p);

  if (lind==p){ // Non-truncated variables
    mean0 = mu; var0 = Sigma; mom20 = var0 + mean0*mean0.t();
  } else {

    if ((lind==0) & (p==1)){ // All variables are truncated: univariate case
      double s11 = as_scalar(sqrt(Sigma));
      double a, b;
      if(lower.is_finite()){ a = as_scalar((lower - mu))/s11; } else { a = -1e40; }
      if(upper.is_finite()){ b = as_scalar((upper - mu))/s11; } else { b = 1e40; }
      double den = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
      double pdfa = R::dnorm4(a,0.0,1.0,0);
      double pdfb = R::dnorm4(b,0.0,1.0,0);
      mean0(0) = as_scalar(mu) + s11*((pdfa - pdfb)/den);
      var0(0,0) = as_scalar(Sigma)*(1.0 + (a*pdfa - b*pdfb)/den - pow((pdfa - pdfb)/den, 2.0));
      mom20(0,0) = var0(0,0) + mean0(0)*mean0(0);
    } else {

      if ((lind==0) & (p>1)){ //All variables are truncated: p-variate case
        arma::mat gen = rtnormal(n,mu,Sigma,lower,upper,burn,thinning);
        mean0 = (mean(gen,0)).t();
        var0 = cov(gen);
        mom20 = var0 + mean0*mean0.t();
      } else {

        if ((lind==(p-1)) & (p>1)){ // One variable is truncated: p-variate case
          double s11 = as_scalar(sqrt(Sigma(ind1,ind1)));
          double a, b;
          if((lower(ind1)).is_finite()){ a = as_scalar((lower(ind1) - mu(ind1))/s11); } else { a = -1e40; }
          if((upper(ind1)).is_finite()){ b = as_scalar((upper(ind1) - mu(ind1))/s11); } else { b = 1e40; }
          double den = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
          double pdfa = R::dnorm4(a,0.0,1.0,0);
          double pdfb = R::dnorm4(b,0.0,1.0,0);
          mean0(ind1) = mu(ind1) + s11*((pdfa - pdfb)/den);
          mean0(ind2) = mu(ind2) + as_scalar((mean0(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var0(ind1,ind1) = Sigma(ind1,ind1)*(1.0 + (a*pdfa - b*pdfb)/den - pow((pdfa - pdfb)/den, 2.0));
          var0(ind2,ind2) = Sigma(ind2,ind2) - ((1.0 - as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1)))/as_scalar(Sigma(ind1,ind1)))*Sigma(ind2,ind1)*Sigma(ind1,ind2);
          var0(ind2,ind1) = as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var0(ind1,ind2) = (var0(ind2,ind1)).t();
          mom20 = var0 + mean0*mean0.t();
        } else {

          if ((lind>0) & (lind<(p-1))){ // The number of truncated variables varies between 2 and p-1
            arma::mat Iden = eye(ind1.size(),ind1.size());
            arma::mat sigInv = Sigma(ind1,ind1).i();
            arma::mat gen = rtnormal(n,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
            mean0(ind1) = (mean(gen,0)).t();
            mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
            var0(ind1,ind1) = cov(gen);
            var0(ind2,ind2) = Sigma(ind2,ind2) - Sigma(ind2,ind1)*sigInv*(Iden - var0(ind1,ind1)*sigInv)*Sigma(ind1,ind2);
            var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
            var0(ind1,ind2) = (var0(ind2,ind1)).t();
            mom20 = var0 + mean0*mean0.t();
          }
        }
      }
    }
  }
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.50*(mom20 + mom20.t());
  output["VarY"] = 0.50*(var0 + var0.t());
  return output;
}


// Compute moments from truncated multivariate Power Exponential distribution
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List PEmoment(arma::vec mu, arma::mat Sigma, double beta, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();
  arma::vec mean0(p); arma::mat mom20(p,p); arma::mat var0(p,p);

  if (lind==p){ // Non-truncated variables
    mean0 = mu;
    var0 = (std::pow(2.0,(1.0/beta))*tgamma((p+2.0)/(2.0*beta))/(p*tgamma(p/(2.0*beta))))*Sigma;
    mom20 = var0 + mean0*mean0.t();

  } else { // There is at least one truncated variable
    arma::mat gen = rtPE(n,beta,mu,Sigma,lower,upper,burn,thinning);
    mean0 = (mean(gen,0)).t();
    var0 = cov(gen);
    mom20 = var0 + mean0*mean0.t();
  }
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.50*(mom20 + mom20.t());
  output["VarY"] = 0.50*(var0 + var0.t());
  return output;
}


// Compute moments from truncated multivariate Pearson VII distribution
// ------------------------------------------------------------------------------------------------------------
List PVIIuni(double mu, double Sigma, double N, double nu, arma::vec lower, arma::vec upper, arma::uword n2){
  double md0 = arma::datum::nan;
  double vr0 = arma::datum::nan;

  if (n2 == 1){
    double s11 = sqrt(Sigma);
    double a = (as_scalar(lower) - mu)/s11;
    double b = (as_scalar(upper) - mu)/s11;
    double nu1 = sqrt((2.0*N - 1.0)/nu);
    double alpha0 = R::pt(nu1*b,2.0*N-1.0,1,0) - R::pt(nu1*a,2.0*N-1.0,1,0);
    double meanx  = 0.0, varx = 0.0;
    if (N == 1.0){
      meanx = sqrt(nu)/(2.0*alpha0*M_PI)*log((nu + b*b)/(nu + a*a));
    } else {
      meanx = sqrt(nu)*tgamma(N)/(2.0*(1-N)*alpha0*sqrt(M_PI)*tgamma(N-0.5))*(pow(1 + b*b/nu, 1.0-N) - pow(1 + a*a/nu, 1-N));
    }
    md0 = mu + s11*meanx; // Mean of Y
    if (N > 1.50){
      double nu2 = sqrt((2.0*N - 3.0)/nu);
      varx = nu*(N-1.0)/(alpha0*(N-1.5))*(R::pt(nu2*b,2.0*N-3.0,1,0) - R::pt(nu2*a,2.0*N-3.0,1,0)) - nu - meanx*meanx;
    } else {
      if (N == 1.50){
        varx = nu/(2.0*alpha0)*(asinh(b/sqrt(nu)) - asinh(a/sqrt(nu))) - nu - meanx*meanx;
      } else {
        varx = nu*tgamma(N)/(alpha0*sqrt(M_PI)*tgamma(N-0.50))*E2integ(1.0-N,a/sqrt(nu),b/sqrt(nu)) - nu - meanx*meanx;
      }
    }
    vr0 = varx*Sigma;
  } else {

    if (N > 1.0){
      double s11 = sqrt(Sigma);
      double nu1 = sqrt((2.0*N-1.0)/nu);
      double a, b;
      if(lower.is_finite()){ a = (as_scalar(lower) - mu)/s11; } else { a = -1e40; }
      if(upper.is_finite()){ b = (as_scalar(upper) - mu)/s11; } else { b = 1e40; }
      double alpha0 = R::pt(nu1*b,2.0*N-1.0,1,0) - R::pt(nu1*a,2.0*N-1.0,1,0);
      double meanx = sqrt(nu)*tgamma(N)/(2.0*(1-N)*alpha0*sqrt(M_PI)*tgamma(N-0.5))*(pow(1 + b*b/nu, 1.0-N) - pow(1 + a*a/nu, 1-N));
      md0 = mu + s11*meanx;
      if (N > 1.50){
        double nu2 = sqrt((2.0*N-3.0)/nu);
        double varx = nu*(N-1.0)/(alpha0*(N-1.5))*(R::pt(nu2*b,2.0*N-3.0,1,0) - R::pt(nu2*a,2.0*N-3.0,1,0)) - nu - meanx*meanx;
        vr0 = varx*Sigma;
      }
    }
  }
  List out;
  out["mean"] = md0; out["var"] = vr0;
  return out;
}

// [[Rcpp::export]]
List PVIImoment(arma::vec mu, arma::mat Sigma, double N, double nu, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::vec mean0(p,fill::zeros);   mean0.replace(0,arma::datum::nan);
  arma::mat var0(p,p,fill::zeros);  var0.replace(0,arma::datum::nan);
  arma::mat mom20(p,p,fill::zeros); mom20.replace(0,arma::datum::nan);
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();

  if (lind == p){ // Non-truncated variables
    if (N > 0.5*(p+1.0)){
      mean0 = mu;
      if (N > 0.5*(p+2.0)){ var0 = nu*Sigma/(2.0*N - p - 2.0); mom20 = var0 + mean0*mean0.t(); }
    }
  } else {

    if ((lind==0) & (p==1)){ // All variables are truncated: univariate case
      arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper));  // Doubly truncated variables
      List momentos = PVIIuni(as_scalar(mu),as_scalar(Sigma),N,nu,lower,upper,ind0.size());
      mean0(0) = as<double>(momentos["mean"]);
      var0(0,0) = as<double>(momentos["var"]);
      mom20(0,0) = var0(0,0) + mean0(0)*mean0(0);
    } else {

      if ((lind==0) & (p>1)){ // All variables are truncated: p-variate case
        arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
        arma::uword n2 = ind0.size();
        arma::uword d = p - n2;
        if (N > (d+2.0)/2.0){ // Exists mean and variance
          arma::mat gen = rtPVII(n,N,nu,mu,Sigma,lower,upper,burn,thinning);
          mean0 = (mean(gen,0)).t();
          var0 = cov(gen); mom20 = var0 + mean0*mean0.t();
        } else {
          if (N > (d+1.0)/2.0){
            if (n2 == 0){ // Exists the mean
              arma::mat gen = rtPVII(n,N,nu,mu,Sigma,lower,upper,burn,thinning);
              mean0 = (mean(gen,0)).t();
            } else { // Exists mean(all), variance(bounded), covariance(all)
              arma::uvec ind3 = unique(join_vert(find_nonfinite(lower),find_nonfinite(upper))); // Non-bounded region
              arma::mat Anan(ind3.size(),ind3.size(),fill::zeros); Anan.replace(0,arma::datum::nan);
              arma::mat gen = rtPVII(n,N,nu,mu,Sigma,lower,upper,burn,thinning);
              mean0 = (mean(gen,0)).t();
              var0 = cov(gen);                   var0(ind3,ind3)  = Anan;
              mom20 = var0 + mean0*mean0.t();    mom20(ind3,ind3) = Anan;
            }
          }
        }
      } else {
        if ((lind==(p-1)) & (p>1)){ // One variable is truncated: p-variate case
          arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
          arma::uword n2 = ind0.size();
          arma::uword d = p - n2;
          if (N > (d+2.0)/2.0){ // Exists mean and variance
            arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper)));   // Truncated variables
            List momentos = PVIIuni(as_scalar(mu(ind1)),as_scalar(Sigma(ind1,ind1)),N-0.5*lind,nu,lower(ind1),upper(ind1),n2);
            arma::vec mx(1);   mx(0) = as<double>(momentos["mean"]);
            arma::mat vx(1,1); vx(0,0) = as<double>(momentos["var"]);
            mean0(ind1) = mx;
            mean0(ind2) = mu(ind2) + as_scalar((mean0(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
            var0(ind1,ind1) = vx;
            double omega21 = (nu + as_scalar((var0(ind1,ind1) + pow(mean0(ind1)-mu(ind1), 2.0))/Sigma(ind1,ind1)))/(2.0*N - lind - 2.0);
            var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - (omega21 - as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1)))/as_scalar(Sigma(ind1,ind1))*Sigma(ind2,ind1)*Sigma(ind1,ind2);
            var0(ind2,ind1) = as_scalar(var0(ind1,ind1)/Sigma(ind1,ind1))*Sigma(ind2,ind1);
            var0(ind1,ind2) = (var0(ind2,ind1)).t();
            mom20 = var0 + mean0*mean0.t();
          } else {
            if (N > (d+1.0)/2.0){
              if (n2 == 0){ // Exists mean
                arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
                List momentos = PVIIuni(as_scalar(mu(ind1)),as_scalar(Sigma(ind1,ind1)),N-0.5*lind,nu,lower(ind1),upper(ind1),n2);
                arma::vec mx(1);   mx(0) = as<double>(momentos["mean"]);
                mean0(ind1) = mx;
                mean0(ind2) = mu(ind2) + as_scalar((mean0(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
              } else { // Exists mean(all), variance(bounded), covariance(all)
                List momentos = PVIIuni(as_scalar(mu(ind0)),as_scalar(Sigma(ind0,ind0)),N-0.5*lind,nu,lower(ind0),upper(ind0),n2);
                arma::vec mx(1);   mx(0) = as<double>(momentos["mean"]);
                arma::mat vx(1,1); vx(0,0) = as<double>(momentos["var"]);
                mean0(ind0) = mx;
                mean0(ind2) = mu(ind2) + as_scalar((mean0(ind0) - mu(ind0))/Sigma(ind0,ind0))*Sigma(ind2,ind0);
                var0(ind0,ind0) = vx;
                var0(ind2,ind0) = as_scalar(var0(ind0,ind0)/Sigma(ind0,ind0))*Sigma(ind2,ind0);
                var0(ind0,ind2) = (var0(ind2,ind0)).t();
                mom20 = var0 + mean0*mean0.t();
              }
            }
          }
        } else { // The number of truncated variables varies between 2 and p-1
          arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
          arma::uvec ind0 = intersect(find_finite(lower),find_finite(upper)); // Doubly truncated variables
          arma::uword n2 = ind0.size();
          arma::uword d = p - n2;

          if (N > (d+2.0)/2.0){ // Exists mean and variance
            arma::mat sigInv = (Sigma(ind1,ind1)).i();
            arma::mat gen = rtPVII(n,(N-0.5*lind),nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
            mean0(ind1) = (mean(gen,0)).t();
            mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
            arma::mat Iden = eye(ind1.size(),ind1.size());
            var0(ind1,ind1) = cov(gen);
            double omega21 = (nu + trace(var0(ind1,ind1)*sigInv) + as_scalar((mean0(ind1) - mu(ind1)).t()*sigInv*(mean0(ind1) - mu(ind1))))/(2.0*N - lind - 2.0);
            var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - Sigma(ind2,ind1)*sigInv*(omega21*Iden - var0(ind1,ind1)*sigInv)*Sigma(ind1,ind2);
            var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
            var0(ind1,ind2) = (var0(ind2,ind1)).t();
            mom20 = var0 + mean0*mean0.t();
          } else {
            if (N > (d+1.0)/2.0){
              if (n2 == 0){ // Exists mean
                arma::mat sigInv = (Sigma(ind1,ind1)).i();
                arma::mat gen = rtPVII(n,(N-0.5*lind),nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
                mean0(ind1) = (mean(gen,0)).t();
                mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
              } else { // Exists mean(all), variance(bounded), covariance(all)
                arma::uvec ind3 = unique(join_vert(find_nonfinite(lower),find_nonfinite(upper))); // Non-bounded region
                arma::mat Anan(ind3.size(),ind3.size(),fill::zeros); Anan.replace(0,arma::datum::nan);
                arma::mat sigInv = (Sigma(ind1,ind1)).i();
                arma::mat gen = rtPVII(n,(N-0.5*lind),nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
                mean0(ind1) = (mean(gen,0)).t();
                mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
                arma::mat Iden = eye(ind1.size(),ind1.size());
                var0(ind1,ind1) = cov(gen);
                var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
                var0(ind1,ind2) = (var0(ind2,ind1)).t();
                mom20 = var0 + mean0*mean0.t();
                var0(ind3,ind3) = Anan; mom20(ind3,ind3) = Anan;
              }
            }
          }
        } // End if
      } // End if
    } // End if
  } // End if
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.50*(mom20 + mom20.t());
  output["VarY"] = 0.50*(var0 + var0.t());
  return output;
}


// Compute moments from truncated multivariate Slash distribution
// ------------------------------------------------------------------------------------------------------------

// w2.1 for the Slash distribution
double wSlash(arma::mat X, arma::vec mu, mat Sinv, double q, double p){
  arma::uword n1 = X.n_rows;
  double cont = 0.0; double delta;
  for (uword i=0;i<n1;i++){
    delta = as_scalar(((X.row(i)).t() - mu).t()*Sinv*((X.row(i)).t() - mu));
    cont += slash_g(delta, q-1.0, p)/slash_g(delta, q, p);
  }
  double omega = cont/n1;
  return omega;
}

// [[Rcpp::export]]
List Slashmoment(arma::vec mu, arma::mat Sigma, double nu, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::vec mean0(p);
  arma::mat var0(p,p,fill::zeros);  var0.replace(0,arma::datum::nan);
  arma::mat mom20(p,p,fill::zeros); mom20.replace(0,arma::datum::nan);
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();

  if (lind==p){ // Non-truncated variables
    mean0 = mu;
    if (nu > 1.0){ var0 = (nu/(nu - 1.0))*Sigma; mom20 = var0 + mean0*mean0.t(); }

  } else {

    if (lind == 0){ // All variables are truncated: p-variate case
      arma::mat gen = rtslash(n,nu,mu,Sigma,lower,upper,burn,thinning);
      mean0 = (mean(gen,0)).t();
      if (nu > 1.0){ var0 = cov(gen); mom20 = var0 + mean0*mean0.t(); }

    } else { // The number of truncated variables varies between 1 and p-1
      arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper)));   // Truncated variables
      arma::mat gen = rtslash(n,nu,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
      arma::mat sigInv = (Sigma(ind1,ind1)).i();
      mean0(ind1) = (mean(gen,0)).t();
      mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
      if (nu > 1.0){
        arma::mat Iden = eye(ind1.size(),ind1.size());
        double omega21 = wSlash(gen, mu(ind1), sigInv, nu, ind1.size());
        var0(ind1,ind1) = cov(gen);
        var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - Sigma(ind2,ind1)*sigInv*(omega21*Iden - var0(ind1,ind1)*sigInv)*Sigma(ind1,ind2);
        var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
        var0(ind1,ind2) = (var0(ind2,ind1)).t();
        mom20 = var0 + mean0*mean0.t();
      }
    }
  }
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.50*(mom20 + mom20.t());
  output["VarY"] = 0.50*(var0 + var0.t());
  return output;
}


// Compute moments from truncated multivariate Contaminated normal distribution
// ------------------------------------------------------------------------------------------------------------

// Compute the density of Multivariate normal distribution
static double const log2pi = std::log(2.0 * M_PI);
arma::vec dmvnorm1(arma::mat const &x,arma::rowvec const &mean,arma::mat const &sigma,bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z      = (x.row(i) - mean) * rooti;
    out(i) = other_terms - 0.5 * arma::dot(z,z);
  }
  if (logd){ return out; }
  return exp(out);
}

// [[Rcpp::export]]
List CNmoment(arma::vec mu, arma::mat Sigma, double nu, double rho, arma::vec lower, arma::vec upper, arma::uword n, int burn, int thinning){

  arma::uword p = mu.size();
  arma::vec mean0(p); arma::mat mom20(p,p); arma::mat var0(p,p);
  arma::uvec ind2 = intersect(find_nonfinite(lower),find_nonfinite(upper));   // Non-truncated variables
  arma::uword lind = ind2.size();

  if (lind==p){ // Non-truncated variables
    mean0 = mu;
    var0 = (nu/rho + 1.0 - nu)*Sigma;
    mom20 = var0 + mean0*mean0.t();
  } else {

    if ((lind==0) & (p==1)){ // All variables are truncated: univariate case
      // Variable 1: X~TN(mu,Sigma/rho)
      double s11 = as_scalar(sqrt(Sigma/rho));
      double a, b;
      if(lower.is_finite()){ a = as_scalar((lower - mu))/s11; } else { a = -1e40; }
      if(upper.is_finite()){ b = as_scalar((upper - mu))/s11; } else { b = 1e40; }
      double den1 = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
      double pdfa = R::dnorm4(a,0.0,1.0,0);
      double pdfb = R::dnorm4(b,0.0,1.0,0);
      double mu1 = as_scalar(mu) + s11*((pdfa - pdfb)/den1);
      double var1 = as_scalar(Sigma)*(1.0 + (a*pdfa - b*pdfb)/den1 - pow((pdfa-pdfb)/den1,2.0))/rho;
      // Variable 2: X~TN(mu,Sigma)
      double s22 = as_scalar(sqrt(Sigma));
      if(lower.is_finite()){ a = as_scalar((lower - mu))/s22; } else { a = -1e40; }
      if(upper.is_finite()){ b = as_scalar((upper - mu))/s22; } else { b = 1e40; }
      double den2 = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
      pdfa = R::dnorm4(a,0.0,1.0,0);
      pdfb = R::dnorm4(b,0.0,1.0,0);
      double mu2 = as_scalar(mu) + s22*((pdfa - pdfb)/den2);
      double var2 = as_scalar(Sigma)*(1.0 + (a*pdfa - b*pdfb)/den2 - pow((pdfa-pdfb)/den2,2.0));
      double pi0 = nu*den1 + (1.0-nu)*den2;
      mean0(0) = (nu*den1*mu1 + (1.0-nu)*den2*mu2)/pi0;
      mom20(0,0) = (nu*den1*(var1 + mu1*mu1) + (1.0-nu)*den2*(var2 + mu2*mu2))/pi0;
      var0(0,0) = mom20(0,0) - mean0(0)*mean0(0);
    } else {

      if ((lind==0) & (p>1)){ // All variables are truncated: p-variate case
        arma::mat gen = rtCN(n,nu,rho,mu,Sigma,lower,upper,burn,thinning);
        mean0 = (mean(gen,0)).t();
        var0 = cov(gen);
        mom20 = var0 + mean0*mean0.t();
      } else {

        if ((lind==(p-1)) & (p>1)){ // One variable is truncated: p-variate case
          arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
          // Variable 1: X~TN(mu,Sigma/rho)
          arma::vec mu1(p);
          arma::mat var1(p,p);
          double s11 = as_scalar(sqrt(Sigma(ind1,ind1)/rho));
          double a, b;
          if((lower(ind1)).is_finite()){ a = as_scalar((lower(ind1) - mu(ind1))/s11); } else { a = -1e40; }
          if((upper(ind1)).is_finite()){ b = as_scalar((upper(ind1) - mu(ind1))/s11); } else { b = 1e40; }
          double den1 = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
          double pdfa = R::dnorm4(a,0.0,1.0,0);
          double pdfb = R::dnorm4(b,0.0,1.0,0);
          mu1(ind1) = mu(ind1) + s11*((pdfa - pdfb)/den1);
          mu1(ind2) = mu(ind2) + as_scalar((mu1(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var1(ind1,ind1) = Sigma(ind1,ind1)*(1.0 + (a*pdfa - b*pdfb)/den1 - pow((pdfa - pdfb)/den1,2.0))/rho;
          var1(ind2,ind2) = Sigma(ind2,ind2)/rho - ((1.0 - as_scalar(rho*var1(ind1,ind1)/Sigma(ind1,ind1)))/as_scalar(Sigma(ind1,ind1)*rho))*Sigma(ind2,ind1)*Sigma(ind1,ind2);
          var1(ind2,ind1) = as_scalar(var1(ind1,ind1)/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var1(ind1,ind2) = (var1(ind2,ind1)).t();
          // Variable 2: X~TN(mu,Sigma)
          arma::vec mu2(p);
          arma::mat var2(p,p);
          s11 = as_scalar(sqrt(Sigma(ind1,ind1)));
          if((lower(ind1)).is_finite()){ a = as_scalar((lower(ind1) - mu(ind1))/s11); } else { a = -1e40; }
          if((upper(ind1)).is_finite()){ b = as_scalar((upper(ind1) - mu(ind1))/s11); } else { b = 1e40; }
          double den2 = R::pnorm5(b,0.0,1.0,1,0) - R::pnorm5(a,0.0,1.0,1,0);
          pdfa = R::dnorm4(a,0.0,1.0,0);
          pdfb = R::dnorm4(b,0.0,1.0,0);
          mu2(ind1) = mu(ind1) + s11*((pdfa - pdfb)/den2);
          mu2(ind2) = mu(ind2) + as_scalar((mu2(ind1) - mu(ind1))/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var2(ind1,ind1) = Sigma(ind1,ind1)*(1.0 + (a*pdfa - b*pdfb)/den2 - pow((pdfa - pdfb)/den2,2.0));
          var2(ind2,ind2) = Sigma(ind2,ind2) - ((1.0 - as_scalar(var2(ind1,ind1)/Sigma(ind1,ind1)))/as_scalar(Sigma(ind1,ind1)))*Sigma(ind2,ind1)*Sigma(ind1,ind2);
          var2(ind2,ind1) = as_scalar(var2(ind1,ind1)/Sigma(ind1,ind1))*Sigma(ind2,ind1);
          var2(ind1,ind2) = (var2(ind2,ind1)).t();
          double pi0 = nu*den1 + (1.0 - nu)*den2;
          mean0 = (den1*nu*mu1 + den2*(1.0 - nu)*mu2)/pi0;
          mom20 = (den1*nu*(var1+mu1*mu1.t()) + den2*(1.0 - nu)*(var2+mu2*mu2.t()))/pi0;
          var0 = mom20 - mean0*mean0.t();

        } else { // The number of truncated variables varies between 2 and p-1
          arma::uvec ind1 = unique(join_vert(find_finite(lower),find_finite(upper))); // Truncated variables
          arma::mat Iden = eye(ind1.size(),ind1.size());
          arma::mat sigInv = (Sigma(ind1,ind1)).i();
          arma::mat gen = rtCN(n,nu,rho,mu(ind1),Sigma(ind1,ind1),lower(ind1),upper(ind1),burn,thinning);
          arma::vec d1 = dmvnorm1(gen,(mu(ind1)).t(),(Sigma(ind1,ind1)/rho),false);
          arma::vec d2 = dmvnorm1(gen,(mu(ind1)).t(),Sigma(ind1,ind1),false);
          double omega1 = mean(nu*d1/(nu*d1 + (1.0-nu)*d2));
          double omega21 = omega1/rho + 1.0 - omega1;
          mean0(ind1) = (mean(gen,0)).t();
          mean0(ind2) = mu(ind2) + Sigma(ind2,ind1)*sigInv*(mean0(ind1) - mu(ind1));
          var0(ind1,ind1) = cov(gen);
          var0(ind2,ind2) = omega21*Sigma(ind2,ind2) - Sigma(ind2,ind1)*sigInv*(omega21*Iden - var0(ind1,ind1)*sigInv)*Sigma(ind1,ind2);
          var0(ind2,ind1) = Sigma(ind2,ind1)*sigInv*var0(ind1,ind1);
          var0(ind1,ind2) = (var0(ind2,ind1)).t();
          mom20 = var0 + mean0*mean0.t();
        }
      }
    }
  }
  List output;
  output["EY"] = mean0;
  output["EYY"] = 0.50*(mom20 + mom20.t());
  output["VarY"] = 0.50*(var0 + var0.t());
  return output;
}
