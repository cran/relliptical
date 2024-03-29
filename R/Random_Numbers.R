#' Sampling Random Numbers from Truncated Multivariate Elliptical Distributions
#'
#' This function generates observations from a truncated multivariate elliptical distribution
#' with location parameter \code{mu}, scale matrix \code{Sigma}, lower and upper
#' truncation points \code{lower} and \code{upper} via Slice Sampling algorithm with Gibbs sampler steps.
#'
#' @param n number of observations to generate. Must be an integer >= 1.
#' @param mu numeric vector of length \eqn{p} representing the location parameter.
#' @param Sigma numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the
#' scale parameter.
#' @param lower vector of lower truncation points of length \eqn{p}.
#' @param upper vector of upper truncation points of length \eqn{p}.
#' @param dist represents the truncated distribution to be used. The values are \code{'Normal'},
#' \code{'t'}, \code{'Laplace'}, \code{'PE'}, \code{'PVII'}, \code{'Slash'}, and \code{'CN'} for the truncated Normal, Student-t,
#' Laplace, Power Exponential, Pearson VII, Slash, and Contaminated Normal distribution, respectively.
#' @param nu additional parameter or vector of parameters depending on the
#' density generating function. See Details.
#' @param expr a character with the density generating function. See Details.
#' @param gFun an R function with the density generating function. See Details.
#' @param ginvFun an R function with the inverse of the density generating function defined in
#' \code{gFun}. See Details.
#' @param burn.in number of samples to be discarded as a burn-in phase.
#' @param thinning factor for reducing the autocorrelation of random points.
#'
#' @details The \code{dist} argument represents the truncated distribution to be used. The values are
#' \code{Normal}, \code{t}, \code{'t'}, \code{PE}, \code{PVII}, \code{Slash}, and \code{CN}, for the
#' truncated Normal, Student-t, Laplace, Power Exponential, Pearson VII, Slash, and Contaminated Normal distribution,
#' respectively.
#'
#' The argument \code{nu} is a parameter or vector of parameters depending on the density generating
#' function (DGF). For the truncated Student-t, Power Exponential, and Slash distribution, \code{nu} is a
#' positive number. For the truncated Pearson VII, \code{nu} is a vector with the first element greater than
#' \eqn{p/2} and the second element a positive number. For the truncated Contaminated Normal distribution,
#' \code{nu} is a vector of length 2 assuming values between 0 and 1.
#'
#' This function also allows generating random numbers from other truncated elliptical distributions not specified
#' in the \code{dist} argument, by supplying the density generating function (DGF) through arguments either
#' \code{expr} or \code{gFun}. The DGF must be a non-negative and strictly decreasing function on \code{(0, Inf)}.
#' The easiest way is to provide the DGF expression to argument \code{expr} as a character. The notation used in
#' \code{expr} needs to be understood by package \code{Ryacas0}, and the environment of \code{R}. For instance,
#' for the DGF \eqn{g(t)=e^{-t}}, the user must provide \code{expr = "exp(1)^(-t)"}. See that the function must depend
#' only on variable \eqn{t}, and any additional parameter must be passed as a fixed value. For this case, when a character
#' expression is provided to \code{expr}, the algorithm tries to compute a closed-form expression for the inverse function
#' of \eqn{g(t)}, however, this is not always possible (a warning message is returned). See example 2.
#'
#' If it was no possible to generate random samples by passing a character expression to \code{expr}, the user may provide
#' a custom \code{R} function to the \code{gFun} argument. By default, its inverse function is approximated numerically,
#' however, the user may also provide its inverse to the \code{ginvFun} argument to gain some computational time.
#' When \code{gFun} is provided, arguments \code{dist} and \code{expr} are ignored.
#'
#' @note The Normal distribution is a particular case of the Power Exponential distribution when \code{nu = 1}.
#' The Student-t distribution with \eqn{\nu} degrees of freedom results from the Pearson VII
#' distribution when \code{nu = } ((\eqn{\nu}+p)/2, \eqn{\nu}).
#'
#' @return It returns a matrix of dimensions \eqn{n}x\eqn{p} with the random points sampled.
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
#' mu  = c(0, 1)
#' Sigma = matrix(c(1,0.70,0.70,3), 2, 2)
#' lower = c(-2, -3)
#' upper = c(3, 3)
#' sample1 = rtelliptical(5e4, mu, Sigma, lower, upper, dist="Normal")
#'
#' # Histogram and density for variable 1
#' ggplot(data.frame(sample1), aes(x=X1)) +
#'    geom_histogram(aes(y=after_stat(density)), colour="black", fill="grey", bins=15) +
#'    geom_density(color="red") + labs(x=bquote(X[1]), y="Density") + theme_bw()
#'
#' # Histogram and density for variable 2
#' ggplot(data.frame(sample1), aes(x=X2)) +
#'    geom_histogram(aes(y=after_stat(density)), colour="black", fill="grey", bins=15) +
#'    geom_density(color="red") + labs(x=bquote(X[2]), y="Density") + theme_bw()
#'
#' \donttest{
#' # Example 2: Sampling from the Truncated Logistic distribution
#'
#' # Function for plotting the sample autocorrelation using ggplot2
#' acf.plot = function(samples){
#'  p = ncol(samples); n = nrow(samples); q1 = qnorm(0.975)/sqrt(n); acf1 = list(p)
#'  for (i in 1:p){
#'    bacfdf = with(acf(samples[,i], plot=FALSE), data.frame(lag, acf))
#'    acf1[[i]] = ggplot(data=bacfdf, aes(x=lag,y=acf)) + geom_hline(aes(yintercept=0)) +
#'      geom_segment(aes(xend=lag, yend=0)) + labs(x="Lag", y="ACF", subtitle=bquote(X[.(i)])) +
#'      geom_hline(yintercept=c(q1,-q1), color="red", linetype="twodash") + theme_bw()
#'  }
#'  return (acf1)
#' }
#'
#' set.seed(5678)
#' mu  = c(0, 0)
#' Sigma = matrix(c(1,0.70,0.70,1), 2, 2)
#' lower = c(-2, -2)
#' upper = c(3, 2)
#' # Sample autocorrelation with no thinning
#' sample2 = rtelliptical(10000, mu, Sigma, lower, upper, expr="exp(1)^(-t)/(1+exp(1)^(-t))^2")
#' grid.arrange(grobs=acf.plot(sample2), top="Logistic distribution with no thinning", nrow=1)
#'
#' # Sample autocorrelation with thinning = 3
#' sample3 = rtelliptical(10000, mu, Sigma, lower, upper, expr="exp(1)^(-t)/(1+exp(1)^(-t))^2",
#'                        thinning=3)
#' grid.arrange(grobs=acf.plot(sample3), top="Logistic distribution with thinning = 3", nrow=1)
#'}
#'
#' # Example 3: Sampling from the Truncated Kotz-type distribution
#' set.seed(5678)
#' mu  = c(0, 0)
#' Sigma = matrix(c(1,-0.5,-0.5,1), 2, 2)
#' lower = c(-2, -2)
#' upper = c(3, 2)
#' sample4 = rtelliptical(2000, mu, Sigma, lower, upper, gFun=function(t){t^(-1/2)*exp(-2*t^(1/4))})
#' f1 = ggplot(data.frame(sample4), aes(x=X1,y=X2)) + geom_point(size=0.50) +
#'      labs(x=expression(X[1]), y=expression(X[2]), subtitle="Kotz(2,1/4,1/2)") + theme_bw()
#' ggMarginal(f1, type="histogram", fill="grey")
#'

#' @references{
#'   \insertRef{fang2018symmetric}{relliptical}
#'
#'   \insertRef{ho2012some}{relliptical}
#'
#'   \insertRef{neal2003slice}{relliptical}
#'
#'   \insertRef{robert2010introducing}{relliptical}
#'
#'   \insertRef{valeriano2023moments}{relliptical}
#' }


rtelliptical = function(n=1e4, mu=rep(0,length(lower)), Sigma=diag(length(lower)), lower,
                        upper=rep(Inf,length(lower)), dist="Normal", nu=NULL, expr=NULL, gFun=NULL,
                        ginvFun=NULL, burn.in=0, thinning=1){

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

  # Validating mu, Sigma, lower and upper dimensions
  if (!(is.vector(lower) | is.matrix(lower)) | any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA")
  if (!(is.vector(upper) | is.matrix(upper)) | any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA")
  if (length(c(lower))!=length(c(upper))) stop ("lower bound must have same dimension than upper")
  if (ncol(as.matrix(lower))>1 | ncol(as.matrix(upper))>1) stop ("lower and upper must have just one column")
  if(!all(lower<upper)) stop ("lower must be smaller than upper")

  if (!all(is.finite(mu)) | !is.numeric(mu) | !(is.vector(mu) | is.matrix(mu))) stop ("mu is not a numeric vector or contains NA")
  if (length(c(lower))!=length(c(mu))) stop ("lower bound must have same dimension than mu")
  if (!all(is.finite(Sigma)) | !is.numeric(Sigma)) stop ("Sigma is not specified or contains NA")
  if (!is.matrix(Sigma)) { Sigma = as.matrix(Sigma) }

  if (length(c(mu)) == 1){
    if (length(c(Sigma)) != 1) stop ("Unconformable dimensions between mu and Sigma")
    if (Sigma <= 0) stop ("Sigma must be a non-negative number")
  } else {
    if (ncol(as.matrix(mu))>1) stop("mu must have just one column")
    if (ncol(as.matrix(Sigma)) != length(c(mu))) stop ("Unconformable dimensions between mu and Sigma")
    if ( !is.symmetric.matrix(Sigma) ){
      stop ("Sigma must be a square symmetrical real positive-definite matrix")
    } else {
      if ( !is.positive.definite(Sigma) ) stop ("Sigma must be a real positive-definite matrix")
    }
  }

  # Generating random numbers
  # ----------------------------------------------------------------------------
  # Validation of gFun
  if (is(gFun,"function")){
    p = length(c(mu))
    quant = qchisq(0.999,df=p)
    g_fun2 = gFun
    if (!is.decreasing(g_fun2, x.bound=c(0,quant), step = 0.1)) stop ("gFun must be a decreasing function")
    if (g_fun2(quant) < 0) stop ("gFun must be a non-negative function")

    # Validation of ginvFun
    if (is(ginvFun,"function")){
      g_inv2 = ginvFun
      eval2 = g_fun2(quant)
      if (!is.decreasing(g_inv2, x.bound=c(0,eval2), step=eval2/10)) stop ("ginvFun must be a decreasing function")
      seque = seq(1,round(quant),length.out=20)
      if (!all(abs(g_inv2(g_fun2(seque)) - c(seque)) < 1e-10)) stop ("ginvFun is not the inverse of gFun")
    } else {
      inverse.fun = function(f){
        function(y){
          uniroot(function(t){f(t) - y}, lower = 0, upper = 1e10, tol=1e-3)[1]$root
        }
      }
      g_inv2 = inverse.fun(gFun)
    }
    out = randomG(n, mu, Sigma, lower, upper, g_fun2, g_inv2, burn.in, thinning)

  } else {

    # Using Ryacas0 for a given function
    if (!is.null(expr)){
      p = length(c(mu))
      quant = qchisq(0.999,df=p)

      # Function y = g(t)
      expr = gsub("\\<x\\>", "t", expr)
      expr = gsub("\\<y\\>", "t", expr)
      expr = gsub("\\<w\\>", "t", expr)
      expr = gsub("\\<z\\>", "t", expr)
      g_fun = yacas(expr)
      g_fun2 = eval(parse(text = paste("function(t){resp = ",g_fun$text,";return(resp)}",sep="")))
      res = try(g_fun2(quant), silent=TRUE)
      if (is(res,"try-error") | res < 0) stop("expr must be a non-negative decreasing function depending only on t")
      if (!is.decreasing(g_fun2, x.bound=c(0,quant), step = 0.01)) stop ("expr must be a decreasing function")

      # Inverse function t = g*(y)
      g_inv = yacas(paste("Solve(", expr, " == y, t)",sep=""))
      mre = length(g_inv$LinAlgForm)
      if (mre == 0) stop("The inverse of expr could not be computed. Try with gFun")
      g_inv1 = strsplit(g_inv$LinAlgForm, "==")[[1]][2]
      g_inv2 = eval(parse(text = paste("function(y){resp =",g_inv1,";return(resp)}",sep="")))
      res2 = try(g_inv2(res), silent=TRUE)
      conta = 2
      while ((is(res2,"try-error") | res2 < 0) & conta <= mre){
        g_inv1 = strsplit(g_inv$LinAlgForm, "==")[[conta]][2]
        g_inv2 = eval(parse(text = paste("function(y){resp =",g_inv1,";return(resp)}",sep="")))
        res2 = try(g_inv2(res), silent=TRUE)
        conta = conta + 1
      }
      if (is(res2,"try-error") | res2 < 0) stop ("The inverse of expr could not be computed. Try with gFun")

      out = randomG(n, mu, Sigma, lower, upper, g_fun2, g_inv2, burn.in, thinning)

    } else {

      if (!is.null(dist)){

        if (dist=="Normal"){
          out = rtnormal(n, mu, Sigma, lower, upper, burn.in, thinning)

        } else {
          if (dist=="t"){
            if (!is.numeric(nu) | length(c(nu))>1){
              stop ("Degrees of freedom (nu) must be a positive number")
            } else {
              if (nu <= 0){
                stop ("Degrees of freedom (nu) must be a positive number")
              } else {
                out = rttrunc(n, nu, mu, Sigma, lower, upper, burn.in, thinning)
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
                  out = rtPE(n, nu, mu, Sigma, lower, upper, burn.in, thinning)
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
                    if (nu[1] <= length(c(mu))/2 | nu[2] <= 0){
                      stop("The first element of nu must be greater than p/2 and the second element must be non-negative")
                    } else {
                      out = rtPVII(n, nu[1], nu[2], mu, Sigma, lower, upper, burn.in, thinning)
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
                      out = rtslash(n, nu, mu, Sigma, lower, upper, burn.in, thinning)
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
                          out = rtCN(n, nu[1], nu[2], mu, Sigma, lower, upper, burn.in, thinning)
                        }
                      }
                    }
                  } else {
                    if (dist=="Laplace"){
                      out = rtPE(n, 0.5, mu, Sigma, lower, upper, burn.in, thinning)
                    } else {
                      stop ("The dist values are `Normal`, `t`, `Laplace`, `PE`, `PVII`, `Slash`, and `CN`")
                    } # End Laplace
                  } # End CN
                } # End Slash
              } # End PVII
            } # End PE
          } # End t
        } # End Normal
      } else {
        stop ("The dist values are `Normal`, `t`, `Laplace`, `PE`, `PVII`, `Slash`, and `CN`, or provide a function through `expr` or `gFun`")
      } # End dist
    } # End expr
  } # End gFun

  return (out)
}
