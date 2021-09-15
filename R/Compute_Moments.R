#' Mean and Variance for Truncated Multivariate Elliptical Distributions
#'
#' This function approximates the mean vector and variance-covariance matrix for some particular truncated elliptical distributions
#' through Monte Carlo integration for the truncated variables and using properties of the conditional expectation for
#' the non-truncated variables. It supports the p-variate Normal (\code{Normal}), Student-t (\code{t}), Power Exponential (\code{PE}),
#' Pearson VII (\code{PVII}), Slash (\code{Slash}), and Contaminated Normal (\code{CN}) distribution.
#'
#' @param mean numeric vector of length \eqn{p} representing the location parameter.
#' @param Sigma numeric positive definite matrix with dimension \eqn{pxp} representing the
#' scale parameter.
#' @param lower vector of lower truncation points of length \eqn{p}.
#' @param upper vector of upper truncation points of length \eqn{p}.
#' @param dist represents the truncated distribution to be used. The values are \code{Normal},
#' \code{t}, \code{PE}, \code{PVII}, \code{Slash} and \code{CN} for the truncated Normal, Student-t,
#' Power Exponential, Pearson VII, Slash and Contaminated Normal distributions, respectively.
#' @param nu additional parameter or vector of parameters depending on the
#' density generating function. See Details.
#' @param n number of Monte Carlo samples to be generated.
#' @param burn.in number of samples to be discarded as a burn-in phase.
#' @param thinning factor for reducing the autocorrelation of random points.
#'
#' @details This function also considers the univariate case. The argument \code{nu} is a parameter
#' or vector of parameters depending on the density generating function (DGF). For the truncated
#' Student-t, Power Exponential, and Slash distribution, \code{nu} is a positive number.
#' For the truncated Pearson VII, \code{nu} is a vector with the first
#' element greater than \eqn{p/2} and the second element a positive number. For the truncated Contaminated Normal
#' distribution, \code{nu} is a vector of length 2 assuming values between 0 and 1.
#'
#' @note The Normal distribution is a particular case of the Power Exponential distribution when \code{nu = 1}.
#' The Student-t distribution with \eqn{\nu} degrees of freedom results from the Pearson VII
#' distribution when \code{nu = } ((\eqn{\nu}+p)/2, \eqn{\nu}).
#'
#' In the Student-t distribution, if \code{nu >= 300}, the Normal case is considered.
#' For Student-t distribution, the algorithm also supports degrees of freedom \code{nu <= 2}.
#' For Pearson VII distribution, the algorithm supports values of \code{m <= (p+2)/2} (first element of \code{nu}).
#'
#' @return It returns a list with three elements:
#' \item{EY}{the mean vector of length \eqn{p}.}
#' \item{EYY}{the second moment matrix of dimensions \eqn{pxp}.}
#' \item{VarY}{the variance-covariance matrix of dimensions \eqn{pxp}.}
#'
#' @author Katherine L. Valeriano, Christian E. Galarza and Larissa A. Matos
#'
#' @seealso \code{\link{rtelliptical}}
#'
#' @examples
#' # Truncated Student-t distribution
#' set.seed(5678)
#' mean = c(0.1, 0.2, 0.3)
#' Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
#'                nrow=length(mean), ncol=length(mean), byrow=TRUE)
#'
#' # Example 1: considering nu = 0.80 and one doubly truncated variable
#' a = c(-0.8, -Inf, -Inf)
#' b = c(0.5, 0.6, Inf)
#' MC11 = mvtelliptical(mean, Sigma, a, b, "t", 0.80)
#'
#' # Example 2: considering nu = 0.80 and two doubly truncated variables
#' a = c(-0.8, -0.70, -Inf)
#' b = c(0.5, 0.6, Inf)
#' MC12 = mvtelliptical(mean, Sigma, a, b, "t", 0.80) # By default n=1e4
#'
#' # Truncated Pearson VII distribution
#' set.seed(9876)
#' MC21 = mvtelliptical(mean, Sigma, a, b, "PVII", c(1.90,0.80), n=1e6) # More precision
#' c(MC12$EY); c(MC21$EY)
#' MC12$VarY;  MC21$VarY
#'
#' # Truncated Normal distribution
#' set.seed(1234)
#' MC31 = mvtelliptical(mean, Sigma, a, b, "Normal", n=1e4)
#' MC32 = mvtelliptical(mean, Sigma, a, b, "Normal", n=1e6) # More precision
#' @references{
#'   \insertRef{fang2018symmetric}{relliptical}
#'
#'   \insertRef{neal2003slice}{relliptical}
#'
#'   \insertRef{robert2010introducing}{relliptical}
#' }
#'
#' @import Ryacas0
#' @import RcppNumerical
#' @importFrom FuzzyNumbers.Ext.2 is.decreasing
#' @importFrom matrixcalc is.positive.definite is.symmetric.matrix
#' @importFrom stats qchisq uniroot
#' @importFrom Rdpack reprompt
#'
#' @export mvtelliptical
#' @export rtelliptical


mvtelliptical = function(mean, Sigma=diag(length(mean)), lower=rep(-Inf,length(mean)),
                         upper=rep(Inf,length(mean)), dist="Normal", nu=NULL, n=1e4,
                         burn.in=0, thinning=3){

  #---------------------------------------------------------------------#
  #                              Validations                            #
  #---------------------------------------------------------------------#

  # Validating mean, Sigma, lower and upper dimensions
  if (is.null(Sigma) | any(is.na(Sigma))) stop ("Sigma is not specified or contains NA")
  if (any(is.na(mean)) | !is.numeric(mean) | !(is.vector(mean) | is.matrix(mean))) stop ("mean is not a numeric vector")
  if (!is.matrix(Sigma)) { Sigma = as.matrix(Sigma) }
  if (!all(c(is.finite(mean)),c(is.finite(Sigma)))) stop ("mean and Sigma must contain only finite values")

  if (length(c(mean)) == 1){
    if (length(c(Sigma)) != 1) stop ("Unconformable dimensions between mean and Sigma")
    if (Sigma <= 0) stop ("Sigma must be a positive real number")
  } else {
    if (ncol(as.matrix(mean))>1) stop("mean must have just one column")
    if (ncol(Sigma) != length(c(mean))) stop ("Unconformable dimensions between mean and Sigma")
    if ( !is.symmetric.matrix(Sigma) ){
      stop ("Sigma must be a square symmetrical real positive-definite matrix")
    } else {
      if ( !is.positive.definite(Sigma) ) stop ("Sigma must be a real positive-definite matrix")
    }
  }
  if (is.null(lower)){
    lower = rep(-Inf,length(mean))
  }else{
    if (length(c(lower))!=length(c(mean))) stop ("lower bound must be numeric and have same dimension than mean")
    if (any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA")
  }
  if (is.null(upper)){
    upper = rep(Inf,length(mean))
  } else {
    if (length(c(upper))!=length(c(mean))) stop ("upper bound must be numeric and have same dimension than mean")
    if (any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA")
  }
  if (all(lower < upper) == FALSE) stop ("lower must be smaller than upper")

  # Validating MCMC parameters
  if (length(c(n))>1 | !is.numeric(n)) stop("n must be an integer")
  if (n%%1!=0 | n<2) stop ("n must be an integer")
  if (length(c(burn.in))>1 | !is.numeric(burn.in)) stop("burn.in must be a non-negative integer")
  if (burn.in%%1!=0 | burn.in<0) stop ("burn.in must be a non-negative integer")
  if (length(c(thinning))>1 | !is.numeric(thinning)) stop("thinning must be a non-negative integer")
  if (thinning%%1!=0 | thinning<1) stop ("thinning must be a non-negative integer")

  if (!is.null(dist)){

    if (dist=="Normal"){
      output = Nmoment(mean, Sigma, lower, upper, n, burn.in, thinning)
    } else {

      if (dist=="t"){
        if (!is.numeric(nu) | length(c(nu))>1){
          stop ("Degrees of freedom (nu) must be a positive number")
        } else {
          if (nu <= 0){
            stop ("Degrees of freedom (nu) must be a positive number")
          } else {
            if (nu < 300){ # t distribution
              output = Tmoment(mean, Sigma, nu, lower, upper, n, burn.in, thinning)
            } else { # Normal distribution
              output = Nmoment(mean, Sigma, lower, upper, n, burn.in, thinning)
            }
          }
        }
      } else {

        if (dist=="PVII"){
          if (!is.numeric(nu)){
            stop ("The vector of parameters nu must be provided for the PVII case")
          } else {
            if (length(c(nu))!=2){
              stop ("The vector of parameters nu must be of length 2")
            } else {
              if (nu[1] <= length(c(mean))/2 | nu[2] <= 0){
                stop("The first element of nu must be greater than p/2 and the second element must be non-negative")
              } else {
                output = PVIImoment(mean, Sigma, nu[1], nu[2], lower, upper, n, burn.in, thinning)
              }
            }
          }
        } else {

          if (dist=="PE"){
            if (!is.numeric(nu) | length(c(nu))>1){
              stop ("Kurtosis parameter (nu) must be a positive number")
            } else {
              if (nu <= 0){
                stop ("Kurtosis parameter (nu) must be a positive number")
              } else {
                if (nu != 1){ # Power exponential distribution
                  output = PEmoment(mean, Sigma, nu, lower, upper, n, burn.in, thinning)
                } else { # Normal distribution
                  output = Nmoment(mean, Sigma, lower, upper, n, burn.in, thinning)
                }
              }
            }
          } else {

            if (dist=="Slash"){
              if (!is.numeric(nu) | length(c(nu))>1){
                stop ("Degrees of freedom (nu) must be non-negative")
              } else {
                if (nu <= 0){
                  stop ("Degrees of freedom (nu) must be non-negative")
                } else {
                  output = Slashmoment(mean, Sigma, nu, lower, upper, n, burn.in, thinning)
                }
              }
            } else {

              if (dist=="CN"){
                if (!is.numeric(nu)){
                  stop ("The vector of parameters nu must be provided for the CN case")
                } else {
                  if (length(c(nu)) != 2){
                    stop ("The vector of parameters nu must be of length 2")
                  } else {
                    if (!all(nu > 0 & nu < 1)){
                      stop ("The values of nu must be between 0 and 1")
                    } else {
                      output = CNmoment(mean, Sigma, nu[1], nu[2], lower, upper, n, burn.in, thinning)
                    }
                  }
                }
              } else {
                stop ("The dist values are Normal, t, PE, PVII, Slash, CN")
              } # End CN
            } # End Slash
          } # End PE
        } # End PVII
      } # End t
    } # End Normal
  } else {
    stop ("The dist values are Normal, t, PE, PVII, Slash, CN")
  }

  return(output)
}
