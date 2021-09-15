#' Sampling Random Numbers from Truncated Multivariate Elliptical Distributions
#'
#' This function generates random numbers from a truncated multivariate elliptical distribution
#' with location parameter equal to \code{mean}, scale matrix \code{Sigma}, lower and upper
#' truncation points \code{lower} and \code{upper} via Slice Sampling algorithm with Gibbs sampler steps.
#'
#' @param n number of observations to generate. Must be an integer >= 1.
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
#' @param expr a character with the density generating function. See Details.
#' @param gFun an R function with the density generating function. See Details.
#' @param ginvFun an R function with the inverse of the density generating function defined in
#' \code{gFun}. See Details.
#' @param burn.in number of samples to be discarded as a burn-in phase.
#' @param thinning factor for reducing the autocorrelation of random points.
#'
#' @details The argument \code{nu} is a parameter or vector of parameters depending on the density generating
#' function (DGF). For the truncated Student-t, Power Exponential, and Slash distribution, \code{nu} is a
#' positive number. For the truncated Pearson VII,
#' \code{nu} is a vector with the first element greater than \eqn{p/2} and the second element a positive number.
#' For the truncated Contaminated Normal distribution, \code{nu} is a vector of length 2 assuming values between 0 and 1.
#'
#' This function also allows us to generate random points from other truncated elliptical distributions by specifying
#' the DGF through arguments \code{expr} or \code{gFun} and making \code{dist} equal to \code{NULL}.
#' The DGF must be non-negative and strictly decreasing on the interval \code{(0, Inf)}.
#' The DGF is given as a character to argument \code{expr}. The notation used in \code{expr} needs to be
#' understood by package \code{Ryacas0} and the environment of R. For example if the DGF is \eqn{g(t)=e^{-t}}, then
#' \code{expr="exp(1)^(-t)"}. In this case, the algorithm tries to compute a closed expression for the inverse function
#' of \eqn{g(t)}, a warning message is returned when it is not possible. Additionally, the function in \code{expr} should
#' depend only on variable \eqn{t}, and any additional parameter must be a fixed value.
#' Argument \code{gFun} can be accessed when \code{expr='NULL'}. It accepts the DGF as an R function.
#' The inverse of the function defined in \code{gFun} can be provided as an R function through \code{ginvFun}.
#' If \code{ginvFun='NULL'}, then the inverse of \code{gFun} is approximated numerically.
#' In the examples, we show how to use \code{expr} and \code{gFun} to draw samples from the bivariate truncated Logistic
#' and Kotz-type distributions. Note that the DGF of Kotz-type distribution is strictly decreasing for values of \eqn{N} between
#' (2-p)/2 and 1, see \insertCite{fang2018symmetric;textual}{relliptical}.
#'
#' @note The Normal distribution is a particular case of the Power Exponential distribution when \code{nu = 1}.
#' The Student-t distribution with \eqn{\nu} degrees of freedom results from the Pearson VII
#' distribution when \code{nu = } ((\eqn{\nu}+p)/2, \eqn{\nu}).
#'
#' @return It returns a matrix of dimensions \eqn{nxp} with the random points sampled.
#'
#' @author Katherine L. Valeriano, Christian E. Galarza and Larissa A. Matos
#'
#' @seealso \code{\link{mvtelliptical}}
#'
#' @examples
#' library(ggplot2)
#' library(ggExtra)
#' library(gridExtra)
#'
#' # Example 1: Sampling from the Truncated Normal distribution
#' set.seed(1234)
#' mean  = c(0, 1)
#' Sigma = matrix(c(1,0.70,0.70,3), 2, 2)
#' lower = c(-2, -3)
#' upper = c(3, 3)
#' sample1 = rtelliptical(5e4, mean, Sigma, lower, upper, dist="Normal")
#'
#' # Histogram and density for variable 1
#' ggplot(data.frame(sample1), aes(x=X1)) +
#'    geom_histogram(aes(y=..density..), colour="black", fill="grey", bins=15) +
#'    geom_density(color="red") + labs(x=bquote(X[1]), y="Density")
#'
#' # Histogram and density for variable 2
#' ggplot(data.frame(sample1), aes(x=X2)) +
#'    geom_histogram(aes(y=..density..), colour="black", fill="grey", bins=15) +
#'    geom_density(color="red") + labs(x=bquote(X[2]), y="Density")
#'
#'
#' # Example 2: Sampling from the Truncated Logistic distribution
#'
#' # Function for plotting the sample autocorrelation using ggplot2
#' acf.plot = function(samples){
#'  p = ncol(samples); n = nrow(samples); q1 = qnorm(0.975)/sqrt(n); acf1 = list(p)
#'  for (i in 1:p){
#'    bacfdf = with(acf(samples[,i], plot=FALSE), data.frame(lag, acf))
#'    acf1[[i]] = ggplot(data=bacfdf, aes(x=lag,y=acf)) + geom_hline(aes(yintercept=0)) +
#'      geom_segment(aes(xend=lag, yend=0)) + labs(x="Lag", y="ACF", subtitle=bquote(X[.(i)])) +
#'      geom_hline(yintercept=c(q1,-q1), color="red", linetype="twodash")
#'  }
#'  return (acf1)
#' }
#'
#' set.seed(5678)
#' mean  = c(0, 0)
#' Sigma = matrix(c(1,0.70,0.70,1), 2, 2)
#' lower = c(-2, -2)
#' upper = c(3, 2)
#' # Sample autocorrelation with no thinning
#' sample2 = rtelliptical(1000, mean, Sigma, lower, upper, dist=NULL,
#'                        expr="exp(1)^(-t)/(1+exp(1)^(-t))^2")
#' grid.arrange(grobs=acf.plot(sample2), top="Logistic distribution with no thinning", nrow=1)
#'
#' # Sample autocorrelation with thinning = 3
#' sample3 = rtelliptical(1000, mean, Sigma, lower, upper, dist=NULL,
#'                        expr="exp(1)^(-t)/(1+exp(1)^(-t))^2", thinning=3)
#' grid.arrange(grobs=acf.plot(sample3), top="Logistic distribution with thinning = 3", nrow=1)
#'
#'
#' # Example 3: Sampling from the Truncated Kotz-type distribution
#' sample4 = rtelliptical(1500, mean, Sigma, lower, upper, dist=NULL, expr=NULL,
#'                        gFun=function(t){ t^(-1/2)*exp(-2*t^(1/4)) })
#' f1 = ggplot(data.frame(sample4), aes(x=X1,y=X2)) + geom_point(size=0.50) +
#'      labs(x=expression(X[1]), y=expression(X[2]), subtitle="Kotz(2,1/4,1/2)")
#' ggMarginal(f1, type="histogram", fill="grey")
#' @references{
#'   \insertRef{fang2018symmetric}{relliptical}
#'
#'   \insertRef{ho2012some}{relliptical}
#'
#'   \insertRef{neal2003slice}{relliptical}
#' }


rtelliptical = function(n=1e4, mean, Sigma=diag(length(mean)), lower=rep(-Inf, length(mean)),
                  upper=rep(Inf, length(mean)), dist="Normal", nu=NULL, expr="exp(1)^(-t)/(1+exp(1)^(-t))^2",
                  gFun=NULL, ginvFun=NULL, burn.in=0, thinning=1){

  #---------------------------------------------------------------------#
  #                              Validations                            #
  #---------------------------------------------------------------------#

  # Validating MCMC parameters
  if (length(c(n))>1 | !is.numeric(n)) stop("n must be an integer")
  if (n%%1!=0 | n<1) stop ("n must be an integer")
  if (length(c(burn.in))>1 | !is.numeric(burn.in)) stop("burn.in must be a non-negative integer")
  if (burn.in%%1!=0 | burn.in<0) stop ("burn.in must be a non-negative integer")
  if (length(c(thinning))>1 | !is.numeric(thinning)) stop("thinning must be a non-negative integer")
  if (thinning%%1!=0 | thinning<1) stop ("thinning must be a non-negative integer")

  # Validating mean, Sigma, lower and upper dimensions
  if (any(is.na(mean)) | !is.numeric(mean) | !(is.vector(mean) | is.matrix(mean))) stop ("mean is not a numeric vector")
  if (is.null(Sigma) | any(is.na(Sigma)) | !is.numeric(Sigma)) stop ("Sigma is not specified or contains NA")
  if (!all(c(is.finite(mean)),c(is.finite(Sigma)))) stop ("mean and Sigma must contain only finite values")
  if (!is.matrix(Sigma)) { Sigma = as.matrix(Sigma) }

  if (length(c(mean)) == 1){
    if (length(c(Sigma)) != 1) stop ("Unconformable dimensions between mean and Sigma")
    if (Sigma <= 0) stop ("Sigma must be a non-negative number")
  } else {
    if (ncol(as.matrix(mean))>1) stop("mean must have just one column")
    if (ncol(as.matrix(Sigma)) != length(c(mean))) stop ("Unconformable dimensions between mean and Sigma")
    if ( !is.symmetric.matrix(Sigma) ){
      stop ("Sigma must be a square symmetrical real positive-definite matrix")
    } else {
      if ( !is.positive.definite(Sigma) ) stop ("Sigma must be a real positive-definite matrix")
    }
  }
  if (is.null(lower)){
    lower = rep(-Inf,length(mean))
  } else {
    if (length(c(lower))!=length(c(mean))) stop ("lower bound must be numeric and have same dimension than mean")
    if (any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA")
  }
  if (is.null(upper)){
    upper = rep(Inf,length(mean))
  } else {
    if (length(c(upper))!=length(c(mean))) stop ("upper bound must be numeric and have same dimension than mean")
    if (any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA")
  }
  if(all(lower<upper) == FALSE) stop ("lower must be smaller than upper")

  # Validating distributions and nu parameter
  if (!is.null(dist)){

    if (dist=="Normal"){
      out = rtnormal(n=n,mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)

    } else {

      if (dist=="t"){
        if (!is.numeric(nu) | length(c(nu))>1){
          stop ("Degrees of freedom (nu) must be a positive number")
        } else {
          if (nu <= 0){
            stop ("Degrees of freedom (nu) must be a positive number")
          } else {
            out = rttrunc(n=n,nu=nu,mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)
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
              out = rtPE(n=n,beta=nu,mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)
            }
          }

        } else {

          if (dist=="PVII"){
            if (!is.numeric(nu)){
              stop ("The vector of parameters nu must be provided for the PVII case")
            } else {
              if (length(c(nu))!= 2){
                stop ("The vector of parameters nu must be of length 2")
              } else {
                if (nu[1] <= length(c(mean))/2 | nu[2] <= 0){
                  stop("The first element of nu must be greater than p/2 and the second element must be non-negative")
                } else {
                  out = rtPVII(n=n,N=nu[1],nu=nu[2],mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)
                }
              }
            }

          } else {

            if (dist=="Slash"){
              if (!is.numeric(nu) | length(c(nu))>1){
                stop ("Degrees of freedom (nu) must be a positive number")
              } else {
                if (nu <= 0){
                  stop ("Degrees of freedom (nu) must be a positive number")
                } else {
                  out = rtslash(n=n,nu=nu,mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)
                }
              }

            } else {
              if (dist == "CN"){
                if (!is.numeric(nu)){
                  stop ("The vector of parameters nu must be provided for the CN case")
                } else {
                  if (length(c(nu))!= 2){
                    stop ("The vector of parameters nu must be of length 2")
                  } else {
                    if (!all(nu > 0 & nu < 1)){
                      stop ("The values of nu must be between 0 and 1")
                    } else {
                      out = rtCN(n=n,nu=nu[1],rho=nu[2],mu=mean,Sigma=Sigma,a=lower,b=upper,burn=burn.in,lag=thinning)
                    }
                  }
                }

              } else {
                stop ("The dist values are Normal, t, PE, PVII, Slash, CN")
              } # End if CN
            } # End if Slash
          } # End if Pearson VII
        } # End if Power exponential
      } # End if t
    } # End if Normal
  } else {

    if ( is.null(expr) & is.null(gFun) ){
      stop ("The dist values are Normal, t, PE, PVII, Slash, CN, or provide a function through expr or gFun")
    } else {

      p = length(c(mean))
      quant = qchisq(0.999,df=p)

      # Using Ryacas0 for a given function
      if (!is.null(expr)){

        # Function y = g(t)
        expr = gsub("\\<x\\>", "t", expr)
        expr = gsub("\\<y\\>", "t", expr)
        expr = gsub("\\<w\\>", "t", expr)
        expr = gsub("\\<z\\>", "t", expr)
        g_fun = yacas(expr)
        g_fun2 = eval(parse(text = paste("function(t){resp = ",g_fun$text,";return(resp)}",sep="")))
        res = try(g_fun2(quant), silent=TRUE)
        if (class(res)=="try-error" | res < 0) stop("expr must be a non-negative decreasing function depending only on t")
        if (!is.decreasing(g_fun2, x.bound=c(0,quant), step = 0.1)) stop ("expr must be a decreasing function")

        # Inverse function t = g*(y)
        g_inv = yacas(paste("Solve(", expr, " == y, t)",sep=""))
        mre = length(g_inv$LinAlgForm)
        if (mre == 0) stop("The inverse of expr could not be computed. Try with gFun")
        g_inv1 = strsplit(g_inv$LinAlgForm, "==")[[1]][2]
        g_inv2 = eval(parse(text = paste("function(y){resp =",g_inv1,";return(resp)}",sep="")))
        res2 = try(g_inv2(res), silent=TRUE)
        conta = 2
        while ((class(res2)=="try-error" | res2 < 0) & conta <= mre){
          g_inv1 = strsplit(g_inv$LinAlgForm, "==")[[conta]][2]
          g_inv2 = eval(parse(text = paste("function(y){resp =",g_inv1,";return(resp)}",sep="")))
          res2 = try(g_inv2(res), silent=TRUE)
          conta = conta + 1
        }
        if (class(res2)=="try-error" | res2 < 0) stop ("The inverse of expr could not be computed. Try with gFun")

      } else {

        # Validation of gFun
        if (class(gFun)=="function"){
          g_fun2 = gFun
          if (!is.decreasing(g_fun2, x.bound=c(0,quant), step = 0.1)) stop ("gFun must be a decreasing function")
          if (g_fun2(quant) < 0) stop ("gFun must be a non-negative function")
        } else {
          stop("gFun is not an R function")
        }
        # Validation of ginvFun
        if (class(ginvFun)=="function"){
          g_inv2 = ginvFun
          if (!is.decreasing(g_inv2, x.bound=c(0,g_fun2(quant)), step = 0.001)) stop ("ginvFun must be a decreasing function")
          if (g_inv2(g_fun2(quant)) < 0) stop ("ginvFun must be a non-negative function")
          seque = seq(1,round(quant),length.out=20)
          if (all(abs(g_inv2(g_fun2(seque)) - c(seque)) < 1e-10) == FALSE) stop ("ginvFun is not the inverse of gFun")
        } else {
          inverse.fun = function(f){
            function(y){
              uniroot(function(t){f(t) - y}, lower = 0, upper = 1e10, tol=1e-3)[1]$root
            }
          }
          g_inv2 = inverse.fun(gFun)
        }
      } # End if

      out = randomG(n=n,mu=mean,Sigma=Sigma,a=lower,b=upper,gFUN=g_fun2,ginvFUN=g_inv2,burn=burn.in,lag=thinning)

    } # End if Ryacas or given function
  } # Else dist == NULL

  return (out)
}
