#' @import stats
#' @importFrom utils glob2rx


# UPDATES 0.2.0.0:
#     -Changing the method of filtering to circular convolution by FFT
#     ---Should be a bit faster, cleans up the code a bit, more accurate when crossing boundary of d=0.5
#     -Added arguments for optimizer
#     -Removed auto-selection of 'include mean', now TRUE by default


m0 <- function(d){
  floor(d+0.5)
}

Phi0 <- function(d){
  d - floor(d+0.5)
}

fft.conv <- function(a, b){
  # Convolve two series, must be the same length
  if (length(a)!=length(b)) stop('Arguments must be vectors of the same length!')
  N <- length(a)
  np2 <- nextn(2*N - 1, 2) # Finds the smallest factor of 2 at least 2N-1, speeds things up
  out <- fft(fft(c(a, numeric(np2 - N))) * fft(c(b, numeric(np2 - N))), inverse = TRUE) / np2
  return(Re(out[1:N]))
}

diff.ar <- function(d, k.max){
  k <- 1:(k.max-1)
  c(1, cumprod((k - d - 1)/k))
}

diff.ma <- function(d, k.max){
  diff.ar(-d, k.max)
}

arma.to.ma <- function(ar=numeric(), ma=numeric(), lag.max){
  # Function to return AR coefficients for ARFIMA model.
  c(1, stats::ARMAtoMA(ar, ma, lag.max))
}

arma.to.ar <- function(ar=numeric(), ma=numeric(), lag.max){
  # Function to return AR coefficients for ARMA model.
  # By symmetry, should be able to just swap and negate inputs to arma.to.ma
  # Trick stolen from astsa::ARMAtoAR function
  arinv <- -1*force(ma)
  mainv <- -1*force(ar)
  c(1, stats::ARMAtoMA(arinv, mainv, lag.max))
}

arma.res <- function(x, ar=numeric(), ma=numeric(), mean.sub=FALSE){
  # Filter a time series with an ARMA model
  if (mean.sub) x <- x-mean(x)
  if (length(ar) + length(ma)>0) {
    out <- fft.conv(arma.to.ar(ar, ma, length(x) - 1), out)
  } else {
    out <- x
  }
  return(out)
}

arfima.res <- function(x, d=0, ar=numeric(0), ma=numeric(0), mean.sub=TRUE){
  if (mean.sub) x <- x - mean(x)

  if (d!=0) {
    m <- m0(d)
    if (m>0){
      x <- diff(c(numeric(m),x), differences = m)
      d <- Phi0(d)
    }
    out <- fft.conv(diff.ar(d, length(x)), x)
  } else {
    out <- x
  }
  if (length(ar) + length(ma)>0) {
    out <- fft.conv(arma.to.ar(ar, ma, length(x) - 1), out)
  }
  return(out)
}

acf.res <- function(y, d=0, ar=numeric(), ma=numeric(),
                    lag.max=floor(sqrt(length(y))), 
                    incl.mean=TRUE){
  # Compute the autocorrelation of the residuals, as defined in Mayoral (2007)
  # incl.mean default behavior assumes there is no deterministic polynomial trend present, see Mayoral (2007)
  # With update, incl.mean should basically always be false. Used to be ifelse(d>=0.5, FALSE, TRUE)

  resids <- arfima.res(y, d, ar, ma, incl.mean)
  resids <- resids - mean(resids) # Not strictly in the algorithm, but if you don't do this you'll get some bad results when incl.mean=FALSE
  denom <- sum(resids^2)
  N <- length(resids)
  
  out <- numeric(lag.max) # prepend zero lag term after the loop, should be equal to denom
  
  for (k in 1:lag.max){
    out[k] <- sum(resids[1:(N-k)]*resids[(1+k):N])
  }
  out <- c(denom, out)/denom
  
  return(out)
}

######

#' Minimum Distance Estimation of ARFIMA Model
#'
#' Fits an ARFIMA(\emph{p},\emph{d},\emph{q}) model to a time series using a minimum distance estimator. For details see Mayoral (2007).
#' @param y Numeric vector of the time series.
#' @param p Autoregressive order.
#' @param q Moving average order.
#' @param d.range Range of allowable values for fractional differencing parameter. Smallest value must be greater than -1.
#' @param start Named vector of length 1 + \code{p} + \code{q} containing initial fit values for the fractional differencing parameter, the AR parameters, and the MA parameters (\emph{e.g.} \code{start = c(d=0.4, ar.1=-0.4, ma.1=0.3, ma.2=0.4)}). If missing, automatically selected.
#' @param lag.max Maximum lag to use when calculating the residual autocorrelations. For details see Mayoral (2007).
#' @param incl.mean Whether or not to include a mean term in the model. The default value of \code{TRUE} is recommended unless the true mean is known and previously subtracted. Mean is returned with standard error, which may be misleading for \code{d>=0.5}.
#' @param verbose Whether to print summary of fit.
#' @param method Method for \code{\link[stats]{optim}}, see \code{help(optim)}.
#' @param control List of additional arguments for \code{\link[stats]{optim}}, see \code{help(optim)}.
#'
#' @return A list containing:\tabular{ll}{
#' \code{pars}\tab A numeric vector of parameter estimates.\cr
#' \tab \cr
#' \code{std.errs} \tab A numeric vector of standard errors on parameters.\cr
#' \tab \cr
#' \code{cov.mat} \tab Parameter covariance matrix (excluding mean).\cr
#' \tab \cr
#' \code{fit.obj} \tab \code{\link[stats]{optim}} fit object.\cr
#' \tab \cr
#' \code{p.val} \tab Ljung-Box \emph{p}-value for fit.\cr
#' \tab \cr
#' \code{residuals} \tab Fit residuals.\cr
#' }
#' @note This method makes no assumptions on the distribution of the innovation series, and the innovation variance does not factor into the covariance matrix of parameter estimates. For this reason, it is not included in the results, but can be estimated from the residuals---see Mayoral (2007).
#' @references Mayoral, L. (2007). Minimum distance estimation of stationary and non-stationary ARFIMA processes. \emph{The Econometrics Journal}, \strong{10}, 124-148. doi: \href{https://doi.org/10.1111/j.1368-423X.2007.00202.x}{10.1111/j.1368-423X.2007.00202.x}
#' @seealso \code{\link{mle.arfima}} for psuedo-maximum likelihood estimation.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- arfima.sim(1000, d=0.6, ar=c(-0.4))
#' fit <- mde.arfima(x, p=1, incl.mean=FALSE, verbose=TRUE) 
#' 
#'
#' ## Fit Summary
#' ## --------------------
#' ## Ljung-Box p-val:  0.276 
#' ##           d     ar.1
#' ## est 0.55031 -0.39261
#' ## err 0.03145  0.03673
#' ##
#' ## Correlation Matrix for ARFIMA Parameters
#' ##           d   ar.1
#' ## d    1.0000 0.6108
#' ## ar.1 0.6108 1.0000
mde.arfima <- function(y, p = 1, q = 0, d.range = c(0, 1), start,
                      lag.max = floor(sqrt(length(y))), incl.mean = TRUE, verbose = FALSE, 
                      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), control = list()){
  
  # Minimum Distance Estimator for ARFIMA models, does not require covariance stationarity
  #
  # Mayoral, L. (2007), Minimum distance estimation of stationary and non-stationary ARFIMA processes.
  # The Econometrics Journal, 10: 124-148. doi:10.1111/j.1368-423X.2007.00202.x
  # 
  # https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1368-423X.2007.00202.x
  
  method <- method[1]
  match.arg(method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"))
  
  # I'm pretty sure incl.d should always be true, it's not clear to me that the estimates are ideal if it's not.
  # Besides, why not just use regular old MLE estimators (e.g. stats::arima) when d=0?
  
  # Adding some more subtle behavior for excluding d without changing function arguments
  incl.d <- TRUE
  if (length(d.range)==1 && d.range[1]==0){
    incl.d <- FALSE
  } else {
    if (length(d.range)!=2) stop("Argument \'d.range\' must be either 0 or a numeric vector of length 2!")
  }
  
  # Checking inputs
  if (any(d.range < -1)) stop("Argument \'d.range\' cannot contain entries smaller than -1!")
  if (!incl.d && p==0 && q==0) stop('Must give non-null model to fit!')
  if ((p %% 1 != 0) || (p %% 1 != 0)) stop('Error: p and q must be positive integers!')
  if ((p < 0) || (q < 0)) stop('Error: p and q must be positive integers!')
  
  # Parameter limits for fit
  ifelse(p>0, ar.upper <- numeric(p)+0.999, ar.upper <- numeric())
  ifelse(p>0, ar.lower <- numeric(p)-0.999, ar.lower <- numeric())
  ifelse(q>0, ma.upper <- numeric(q)+0.999, ma.upper <- numeric())
  ifelse(q>0, ma.lower <- numeric(q)-0.999, ma.lower <- numeric())

  d.upper <- max(d.range)
  d.lower <- min(d.range)
  
  upper <- c(d.upper, ar.upper, ma.upper)
  lower <- c(d.lower, ar.lower, ma.lower)
  
  if (!incl.d){
    upper <- upper[-1]
    lower <- lower[-1]
  }
  
  # To avoid initializing parameters at zero. Can result in slight variation in location of fit convergence.
  # Includes name scheme to make things easier later.
  if (missing(start)){
    start <- (upper + lower)/2 + runif(length(upper), -.5, 0.5)
    name.vec <- character()
    if (incl.d) name.vec <- c('d')
    if (p > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ar.',i), 1:p))
    if (q > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ma.',i), 1:q))
    names(start) <- name.vec
  }
  # This is the function to be minimized, sum of the squares of the residual autocorrelations
  C <- function(pars){
    if (any(pars > upper) || any(pars < lower)) return(Inf)
    
    # If pars are named, this should do it.
    d.par <- pars[grep(glob2rx("d"), names(pars), value=TRUE)]
    ar.pars <- pars[grep(glob2rx("ar*"), names(pars), value=TRUE)]
    ma.pars <- pars[grep(glob2rx("ma*"), names(pars), value=TRUE)]
    
    if (length(d.par)==0) d.par <- c(d=0)
    
    # Currently, this overrides the default behavior of acf.res for incl.mean
    # This is to avoid fitting different functions for if d lands on zero in the optimization
    sum((acf.res(y, d.par, ar.pars, ma.pars, lag.max, incl.mean)^2))
  }

  fit <- optim(start, fn=C, method=method, control=control)
  
  # Picking apart fit$par to get estimates to use in finding standard error.
  
  if (incl.d){
    d.par <- fit$par['d']
  } else {
    d.par <- 0
  }
  ar.pars <- fit$par[grep(glob2rx("ar*"), names(fit$par), value=TRUE)]
  ma.pars <- fit$par[grep(glob2rx("ma*"), names(fit$par), value=TRUE)]


  # Goodness-of-fit check.
  # unname just in case, to avoid names on parameters carrying through
  # Max lag value suggestion is from Mayoral (2007)
  resids <- unname(acf.res(y, d=d.par, ar=ar.pars, ma=ma.pars, incl.mean = incl.mean, lag.max=floor(length(y)^0.25)))
  k.vec <- 1:floor(length(y)^0.25)
  sc.resids <- resids[-1]^2/(length(y) - k.vec)
  TS <- length(y)*(length(y)+2)*sum(sc.resids)
  dof <- floor(length(y)^0.25) - p - q
  p.val <- pchisq(TS, dof)
  
  # Computing mean value and standard error
  if (incl.mean){
    if (incl.d){
      # if (m0(d.par) > 0){ ######################################## This depends on a subtle definition of 'mean', ignore it.
      #   mu <- mean(diff(y, differences = m0(d.par)))
      #   se.mu <- sd(diff(y, differences = m0(d.par)))/length(diff(y, differences = m0(d.par)))
      # } else {
      #   mu <- mean(y)
      #   se.mu <- sd(y)/length(y)
      # }
      mu <- mean(y)
      se.mu <- sd(y)/sqrt(length(y))
    } else {
      mu <- mean(y)
      se.mu <- sd(y)/sqrt(length(y))
    }
  }
  
  # Computing standard errors requires polynomial expansion of inverse AR/MA polynomials
  # To get omega weights, want MA coefficients for AR part of model
  # To get psi weights, want AR coefficients for MA part of model
  J <- matrix(0, nrow=lag.max, ncol=ifelse(incl.d,1+q+p, q+p)) # Columns correspond to non-mean model parameters
  
  # Construct entries for matrix J used to estimate standard errors, see Mayoral (2007)
  if (incl.d){
    J[,1] <- -1/(1:lag.max)
    if (p > 0){
      for (i in 1:p){
        J[i:lag.max, i+1] <- arma.to.ma(ar=ar.pars, lag.max = lag.max - i)
      }
      if (q > 0){
        for (i in 1:q){
          J[i:lag.max, p+1+i] <- arma.to.ar(ma=ma.pars, lag.max = lag.max - i)
        }
      }
    } else {
      if (q > 0){
        for (i in 1:q){
          J[i:lag.max, 1+i] <- arma.to.ar(ma=ma.pars, lag.max = lag.max - i)
        }
      }
    }
  } else {
    if (p > 0){
      for (i in 1:p){
        J[i:lag.max, i] <- arma.to.ma(ar=ar.pars, lag.max = lag.max - i)
      }
      if (q > 0){
        for (i in 1:q){
          J[i:lag.max, p+i] <- arma.to.ar(ma=ma.pars, lag.max = lag.max - i)
        }
      }
    } else {
      if (q > 0){
        for (i in 1:q){
          J[i:lag.max, i] <- arma.to.ar(ma=ma.pars, lag.max = lag.max - i)
        }
      }
    }
  }
  
  Xi <- t(J) %*% J
  cov.mat <- solve(Xi)/length(y)
  
  std.errs <- sqrt(diag(cov.mat))
  
  # Picking apart covariances to pretty print standard errors
  # Sorry this is ugly, I'm too lazy to see if I can make it cleaner using parameter names
  if (incl.d){
    d.err <- std.errs[1]
    if (p > 0){
      ar.errs <- std.errs[2:(p+1)]
      if (q > 0) ma.errs <- std.errs[(p+2):length(std.errs)]
    } else {
      if (q > 0) ma.errs <- std.errs[2:length(std.errs)]
    }
  } else {
    if (p > 0){
      ar.errs <- std.errs[1:p]
      if (q > 0) ma.errs <- std.errs[(p+1):length(std.errs)]
    } else {
      if (q > 0) ma.errs <- std.errs[1:length(std.errs)]
    }
  }
  
  pars <- fit$par
  name.vec <- character()
  if (incl.d) name.vec <- c('d')
  if (p > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ar.',i), 1:p))
  if (q > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ma.',i), 1:q))
  rownames(cov.mat) <- name.vec
  colnames(cov.mat) <- name.vec
  if (incl.mean){
    name.vec <- c('mu', name.vec)
    pars <- c(mu, pars)
    std.errs <- c(se.mu, std.errs)
  }
  names(pars) <- name.vec
  names(std.errs) <- name.vec
  
  
  res <- arfima.res(y, d=d.par, ar=ar.pars, ma=ma.pars, mean.sub = incl.mean)
  
  out <- list(pars=pars, std.errs=std.errs, cov.mat=cov.mat, fit.obj=fit, p.val=p.val, residuals=res)
  
  if (verbose){
    cat('\n\nFit Summary')
    cat('\n--------------------')
    cat('\nLjung-Box p-val: ', signif(p.val, 3),'\n')
    estim.frame <- data.frame(tmp=c(0,0))
    if (incl.mean){
      estim.frame <- cbind(estim.frame, c(mu, se.mu))
    }
    if (incl.d){
      estim.frame <- cbind(estim.frame, c(d.par, d.err))
    }

    if (p > 0){
      ar.mat <- t(matrix(c(ar.pars, ar.errs), ncol=2))
      ar.frame <- as.data.frame(ar.mat)
      colnames(ar.frame) <- mapply(function(i) paste0('ar.',i), 1:p)
      estim.frame <- cbind(estim.frame, ar.frame)
    }
    if (q > 0){
      ma.mat <- t(matrix(c(ma.pars, ma.errs), ncol=2))
      ma.frame <- as.data.frame(ma.mat)
      colnames(ma.frame) <- mapply(function(i) paste0('ma.',i), 1:q)
      estim.frame <- cbind(estim.frame, ma.frame)
    }
    estim.frame <- data.frame(estim.frame[,-1]) # data.frame wrapping is to fix bug that happens when p=q=0 and incl.mean=FALSE, otherwise returns a vector and screws up naming step
    colnames(estim.frame) <- name.vec
    rownames(estim.frame) <- c('est','err')
    print(format(estim.frame, digits=4))
    if (ifelse(incl.d, 1,0) + p + q > 1){
      cat('\nCorrelation Matrix for ARFIMA Parameters\n')
      print(format(as.data.frame(cov2cor(cov.mat)), digits=4))}
  }
  
  return(out)
}



#' Pseudo-Maximum Likelihood Estimation of ARFIMA Model
#'
#' Fits an ARFIMA(\emph{p},\emph{d},\emph{q}) model to a time series using a pseudo-maximum likelihood estimator. For details see Beran (1995).
#' @param y Numeric vector of the time series.
#' @param p Autoregressive order.
#' @param q Moving average order.
#' @param d.range Range of allowable values for fractional differencing parameter. Smallest value must be greater than -1.
#' @param start Named vector of length 1 + \code{p} + \code{q} containing initial fit values for the fractional differencing parameter, the AR parameters, and the MA parameters (\emph{e.g.} \code{start = c(d=0.4, ar.1=-0.4, ma.1=0.3, ma.2=0.4)}). If missing, automatically selected.
#' @param incl.mean Whether or not to include a mean term in the model. The default value of \code{TRUE} is recommended unless the true mean is known and previously subtracted. Mean is returned with standard error, which may be misleading for \code{d>=0.5}.
#' @param verbose Whether to print summary of fit.
#' @param method Method for \code{\link[stats]{optim}}, see \code{help(optim)}.
#' @param control List of additional arguments for \code{\link[stats]{optim}}, see \code{help(optim)}.
#'
#' @return A list containing:\tabular{ll}{
#' \code{pars}\tab A numeric vector of parameter estimates.\cr
#' \tab \cr
#' \code{std.errs} \tab A numeric vector of standard errors on parameters.\cr
#' \tab \cr
#' \code{cov.mat} \tab Parameter covariance matrix (excluding mean).\cr
#' \tab \cr
#' \code{fit.obj} \tab \code{\link[stats]{optim}} fit object.\cr
#' \tab \cr
#' \code{p.val} \tab Ljung-Box \emph{p}-value for fit.\cr
#' \tab \cr
#' \code{residuals} \tab Fit residuals.\cr
#' }
#' @references Beran, J. (1995). Maximum Likelihood Estimation of the Differencing Parameter for Short and Long Memory Autoregressive Integrated Moving Average Models. \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \strong{57}, No. 4, 659-672. doi: \href{https://doi.org/10.1111/j.2517-6161.1995.tb02054.x}{10.1111/j.2517-6161.1995.tb02054.x}
#' @seealso \code{\link{mde.arfima}} for minimum distance estimation.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- arfima.sim(1000, d=0.6, ar=c(-0.4))
#' fit <- mle.arfima(x, p=1, incl.mean=FALSE, verbose=TRUE) 
#' 
#' ## Fit Summary
#' ## --------------------
#' ## Ljung-Box p-val:  0.263 
#' ##        sig2       d     ar.1
#' ## est 1.09351 0.53933 -0.37582
#' ## err 0.05343 0.04442  0.05538
#' ##
#' ## Correlation Matrix for ARFIMA Parameters
#' ##         sig2       d    ar.1
#' ## sig2  1.0000 -0.3349  0.4027
#' ## d    -0.3349  1.0000 -0.8318
#' ## ar.1  0.4027 -0.8318  1.0000
mle.arfima <- function(y, p = 1, q = 0, d.range = c(0, 1), start, incl.mean = TRUE, verbose = FALSE, 
                       method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), control = list()){
  # Psuedo-maximum likelihood Estimator for ARFIMA models, does not require covariance stationarity
  #
  # Beran, J. (1995), Maximum Likelihood Estimation of the Differencing Parameter for Invertible Short and Long Memory Autoregressive Integrated Moving Average Models.
  # Journal of the Royal Statistical Society: Series B (Methodological), 57: 659-672. doi:10.1111/j.2517-6161.1995.tb02054.x
  #
  # https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1995.tb02054.x
  
  method <- method[1]
  match.arg(method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"))
  
  # Basic idea is simple: filter the series with an AR filter to get residuals, minimuze the sum of squares
  # Mostly the same as MDE estimator, most code copied.

  # I'm pretty sure incl.d should always be true, it's not clear to me that the estimates are ideal if it's not.
  # Besides, why not just use regular old MLE estimators (e.g. stats::arima) when d=0?
  
  # Adding some more subtle behavior for excluding d without changing function arguments
  incl.d <- TRUE
  if (length(d.range)==1 && d.range[1]==0){
    incl.d <- FALSE
  } else {
    if (length(d.range)!=2) stop("Argument \'d.range\' must be either 0 or a numeric vector of length 2!")
  }

  # Checking inputs
  if (any(d.range < -1)) stop("Argument \'d.range\' cannot contain entries smaller than -1!")
  if (!incl.d && p==0 && q==0) stop('Must give non-null model to fit!')
  if ((p %% 1 != 0) || (p %% 1 != 0)) stop('Error: p and q must be positive integers!')
  if ((p < 0) || (q < 0)) stop('Error: p and q must be positive integers!')
  
  # Parameter limits for fit
  ifelse(p>0, ar.upper <- numeric(p)+0.999, ar.upper <- numeric())
  ifelse(p>0, ar.lower <- numeric(p)-0.999, ar.lower <- numeric())
  ifelse(q>0, ma.upper <- numeric(q)+0.999, ma.upper <- numeric())
  ifelse(q>0, ma.lower <- numeric(q)-0.999, ma.lower <- numeric())

  d.upper <- max(d.range)
  d.lower <- min(d.range)
  
  upper <- c(d.upper, ar.upper, ma.upper)
  lower <- c(d.lower, ar.lower, ma.lower)

  if (!incl.d){
    upper <- upper[-1]
    lower <- lower[-1]
  }
  
  # To avoid initializing parameters at zero. Can result in slight variation in location of fit convergence.
  # Includes name scheme to make things easier later.
  if (missing(start)){
    start <- (upper + lower)/2 + runif(length(upper), -.5, 0.5)
    name.vec <- character()
    if (incl.d) name.vec <- c('d')
    if (p > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ar.',i), 1:p))
    if (q > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ma.',i), 1:q))
    names(start) <- name.vec
  }
  
  
  C <- function(pars){
    if (any(pars > upper) || any(pars < lower)) return(Inf)
    
    # If pars are named, this should do it.
    d.par <- pars[grep(glob2rx("d"), names(pars), value=TRUE)]
    ar.pars <- pars[grep(glob2rx("ar*"), names(pars), value=TRUE)]
    ma.pars <- pars[grep(glob2rx("ma*"), names(pars), value=TRUE)]
    
    if (length(d.par)==0) d.par <- c(d=0)
    
    sum((arfima.res(y, d.par, ar.pars, ma.pars, incl.mean)[-1]^2))
  }
  
  fit <- optim(start, fn=C, method=method, control=control)
  
  # Picking apart fit$par to get estimates to use in finding standard error.
  if (incl.d){
    d.par <- fit$par['d']
  } else {
    d.par <- 0
  }
  ar.pars <- fit$par[grep(glob2rx("ar*"), names(fit$par), value=TRUE)]
  ma.pars <- fit$par[grep(glob2rx("ma*"), names(fit$par), value=TRUE)]
  
  # Goodness-of-fit check.
  # unname just in case, to avoid names on parameters carrying through
  # Max lag value suggestion is from Mayoral (2007)
  resids <- unname(acf.res(y, d=d.par, ar=ar.pars, ma=ma.pars, incl.mean = incl.mean, lag.max=floor(length(y)^0.25)))
  k.vec <- 1:floor(length(y)^0.25)
  sc.resids <- resids[-1]^2/(length(y) - k.vec)
  TS <- length(y)*(length(y)+2)*sum(sc.resids)
  dof <- floor(length(y)^0.25) - p - q
  p.val <- pchisq(TS, dof)
  
  # Computing mean value and standard error
  if (incl.mean){
    if (incl.d){
      # if (m0(d.par) > 0){ ######################################## This depends on a subtle definition of 'mean', ignore it.
      #   mu <- mean(diff(y, differences = m0(d.par)))
      #   se.mu <- sd(diff(y, differences = m0(d.par)))/length(diff(y, differences = m0(d.par)))
      # } else {
      #   mu <- mean(y)
      #   se.mu <- sd(y)/length(y)
      # }
      mu <- mean(y)
      se.mu <- sd(y)/sqrt(length(y))
    } else {
      mu <- mean(y)
      se.mu <- sd(y)/sqrt(length(y))
    }
  }
  
  # Computing innovation variance
  sig2 <- C(fit$par)/(length(y)-1-m0(d.par))
  
  
  # Alrighty, then. V=2*solve(D) gives the covariance of the parameter estimates
  # D_ij given by (2*pi)^(-1)*integrate(dlogspec(x; theta*), -pi, pi)
  # theta* is the parameter vector with only the fractional part of d: (sigma2, Phi0(d), ar.pars, ma.pars)
  # dlogspec is the product of the partial derivatives of the log spectrum with respect to theta_i and theta_j
  
  # AR/MA polynomials and their derivatives
  ar.poly <- function(x, phi){
    p <- length(phi)
    1 - sum(mapply(function(i) phi[i]*exp(1i*x*i), 1:p))
  }
  dar.poly <- function(x, phi, l){
    -2*Re(ar.poly(x, phi)*(-1*exp(l*1i*x)))/Mod(ar.poly(x, phi))
  }
  
  ma.poly <- function(x, theta){
    q <- length(theta)
    1 + sum(mapply(function(i) theta[i]*exp(1i*x*i), 1:q))
  }
  dma.poly <- function(x, theta, l){
    2*Re(ma.poly(x, theta)*(1*exp(l*1i*x)))/Mod(ma.poly(x, theta))
  }


  dlogspec <- function(x, sigma2, pars, i, j){
    if (incl.d){
      d.par <- pars[1]
      if (p > 0){
        ar.pars <- pars[2:(p+1)]
        if (q > 0){
          ma.pars <- pars[(p+2):length(pars)]
        } else {
          ma.pars <- numeric()
        }
      } else {
        ar.pars <- numeric()
        if (q > 0){
          ma.pars <- pars[2:length(pars)]
        } else {
          ma.pars <- numeric()
        }
      }
    } else {
      d.par <- 0
      if (p > 0){
        ar.pars <- pars[1:p]
        if (q > 0){
          ma.pars <- pars[(p+1):length(pars)]
        } else {
          ma.pars <- numeric()
        }
      } else {
        ar.pars <- numeric()
        if (q > 0){
          ma.pars <- pars[1:length(pars)]
        } else {
          ma.pars <- numeric() # This should not happen, but just in case.
        }
      }
    }
    
    d.par <- m0(d.par)
    
    # First do the i part
    if (i==1){
      i.part <- 1/sigma2
    } else if (i==2) {
      i.part <- -2*log(Mod(1-exp(1i*x)))
    } else if (i <= 2 + length(ar.pars)){
      # This means the derivative should be of the ar part of the log spectrum
      i.part <- dar.poly(x, ar.pars, i-2)
    } else {
      i.part <- dma.poly(x, ma.pars, i-2-length(ar.pars))
    }
    
    # Now do the j part
    if (j==1){
      j.part <- 1/sigma2
    } else if (j==2) {
      j.part <- -2*log(Mod(1-exp(1i*x)))
    } else if (j <= 2 + length(ar.pars)){
      # This means the derivative should be of the ar part of the log spectrum
      j.part <- dar.poly(x, ar.pars, j-2)
    } else {
      j.part <- dma.poly(x, ma.pars, j-2-length(ar.pars))
    }

    return(i.part*j.part)
    
  }
  vec.dlogspec <- Vectorize(dlogspec, 'x')
  # This is a symmetric matrix, below only returns diagonal and lower triangular part
  D.mat <- matrix(0, nrow=1+length(fit$par), ncol=1+length(fit$par))

  for (i1 in 1:(1+length(fit$par))){
    for (j1 in 1:i1){
      if (i1==2 && j1==2){
        D.mat[i1,j1] <- (pi^2)/3
      } else {
        # These are symmetric integrals that diverge at the origin, changed bounds and multiplied by 2.
        D.mat[i1,j1] <- integrate(vec.dlogspec, lower=0, upper=pi, sigma2=sig2, pars=fit$par, i=i1, j=j1)$value/pi
      }
    }
  }
  # Constructing full matrix
  D.mat <- D.mat + t(D.mat - diag(diag(D.mat)))
  var.mat <- 2*solve(D.mat)/length(y)
  
  std.errs <- sqrt(diag(var.mat))
  
  
  # Picking apart covariances to pretty print standard errors
  sig2.err <- std.errs[1]
  if (incl.d){
    d.err <- std.errs[2]
    if (p > 0){
      ar.errs <- std.errs[3:(p+2)]
      if (q > 0) ma.errs <- std.errs[(p+3):length(std.errs)]
    } else {
      if (q > 0) ma.errs <- std.errs[3:length(std.errs)]
    }
  } else {
    if (p > 0){
      ar.errs <- std.errs[2:p]
      if (q > 0) ma.errs <- std.errs[(p+2):length(std.errs)]
    } else {
      if (q > 0) ma.errs <- std.errs[2:length(std.errs)]
    }
  }
  
  pars <- c(sig2, fit$par)
  name.vec <- c('sig2')
  if (incl.d) name.vec <- c(name.vec, 'd')
  if (p > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ar.',i), 1:p))
  if (q > 0) name.vec <- c(name.vec, mapply(function(i) paste0('ma.',i), 1:q))
  rownames(var.mat) <- name.vec
  colnames(var.mat) <- name.vec
  if (incl.mean){
    name.vec <- c('mu', name.vec)
    pars <- c(mu, pars)
    std.errs <- c(se.mu, std.errs)
  }
  names(pars) <- name.vec
  names(std.errs) <- name.vec
  # 
  # # Goodness-of-fit based on correlation of residuals
  # resids <- unname(acf.res(y, d=d.par, ar=ar.pars, ma=ma.pars, incl.mean = incl.mean, lag.max=floor(length(y)^0.25)))
  # k.vec <- 1:floor(length(y)^0.25)
  # sc.resids <- resids[-1]^2/(length(y) - k.vec)
  # TS <- length(y)*(length(y)+2)*sum(sc.resids)
  # dof <- floor(length(y)^0.25) - p - q
  # p.val <- pchisq(TS, dof)
  # 
  res <- arfima.res(y, d=d.par, ar=ar.pars, ma=ma.pars, mean.sub = incl.mean)
  
  out <- list(pars=pars, std.errs=std.errs, cov.mat=var.mat, fit.obj=fit, p.val=p.val, residuals=res)
  
  if (verbose){
    cat('\n\nFit Summary')
    cat('\n--------------------')
    cat('\nLjung-Box p-val: ', signif(p.val, 3),'\n')
    estim.frame <- data.frame(sig2=c(sig2, sig2.err))
    if (incl.mean){
      estim.frame <- cbind(c(mu, se.mu), estim.frame)
    }
    if (incl.d){
      estim.frame <- cbind(estim.frame, c(d.par, d.err))
    }
    
    if (p > 0){
      ar.mat <- t(matrix(c(ar.pars, ar.errs), ncol=2))
      ar.frame <- as.data.frame(ar.mat)
      colnames(ar.frame) <- mapply(function(i) paste0('ar.',i), 1:p)
      estim.frame <- cbind(estim.frame, ar.frame)
    }
    if (q > 0){
      ma.mat <- t(matrix(c(ma.pars, ma.errs), ncol=2))
      ma.frame <- as.data.frame(ma.mat)
      colnames(ma.frame) <- mapply(function(i) paste0('ma.',i), 1:q)
      estim.frame <- cbind(estim.frame, ma.frame)
    }
    
    colnames(estim.frame) <- name.vec
    rownames(estim.frame) <- c('est','err')
    print(format(estim.frame, digits=4))
    if (ifelse(incl.d, 1,0) + p + q > 0){
      cat('\nCorrelation Matrix for ARFIMA Parameters\n')
      print(format(as.data.frame(cov2cor(var.mat)), digits=4))}
  }
  
  return(out)
  
}

#' Simulate ARFIMA Process
#'
#' Simulates a series under the given ARFIMA model by applying an MA filter to a series of innovations.
#' 
#' The model is defined by values for the AR and MA parameters (\eqn{\phi} and \eqn{\theta}, respectively), along with the fractional differencing parameter \emph{d}. When \eqn{d\geq 0.5}, then the integer part is taken as \eqn{m=\lfloor d+0.5\rfloor}, and the remainder (between -0.5 and 0.5) stored as \emph{d}. For \eqn{m=0}, the model is:
#' \deqn{\left(1 - \sum_{i=1}^p \phi_i B^i\right)\left(1 - B\right)^d (y_t - \mu)=\left(1 + \sum_{i=1}^q \theta_i B^i\right) \epsilon_t}
#' where \emph{B} is the backshift operator (\eqn{B y_t = y_{t-1}}) and \eqn{\epsilon_t} is the innovation series. When \eqn{m > 0}, the model is defined by:
#' \deqn{y_t = (1 - B)^{-m}x_t}
#' \deqn{\left(1 - \sum_{i=1}^p \phi_i B^i\right)(1 - B)^d (x_t - \mu)=\left(1 + \sum_{i=1}^q \theta_i B^i\right) \epsilon_t}
#' When \code{stat.int = FALSE}, the differencing filter applied to the innovations is not split into parts, and the series model follows the first equation regardless of the value of \emph{d}. This means that \eqn{\mu} is added to the series after filtering and before any truncation. When \code{stat.int = TRUE}, \eqn{x_t - \mu} is generated from filtered residuals, \eqn{\mu} is added, and the result is cumulatively summed \emph{m} times. For non-zero mean and \eqn{m>0}, this will yield a polynomial trend in the resulting data.
#' 
#' Note that the burn-in length may affect the distribution of the sample mean, variance, and autocovariance. Consider this when generating ensembles of simulated data
#' @param n Desired series length.
#' @param d Fractional differencing parameter.
#' @param ar Vector of autoregressive parameters.
#' @param ma Vector of moving average parameters, following the same sign convention as \code{\link[stats]{arima}}.
#' @param mu Mean of process. By default, added after integer integration but before burn-in truncation (see \code{stat.int}).
#' @param sig2 Innovation variance if innovations not provided.
#' @param stat.int Controls integration for non-stationary values of \code{d} (\emph{i.e.} \code{d>=0.5}). If \code{TRUE}, \code{d} split into integer part and stationary part, which will result in a trend when \code{d>=0.5} and \code{mu!=0}.
#' @param n.burn Number of burn-in steps. If not given, chosen based off presence of long memory (\emph{i.e.} \code{d>0}).
#' @param innov Series of innovations. Drawn from normal distribution if not given.
#' @param exact.innov Whether to force the exact innovation series to be used. If \code{FALSE}, innovations will be prepended with resampled points as needed to match \code{n+n.burn}.
#'
#' @return A numeric vector of length n.
#' @export
#'
#' @examples
#' ## Generate ARFIMA(1,d,0) series with Gaussian innovations
#' x <- arfima.sim(1000, d=0.6, ar=c(-0.4)) 
#' 
#' ## Generate ARFIMA(1,d,0) series with uniform innovations.
#' innov.series <- runif(1000, -1, 1)
#' x <- arfima.sim(1000, d=0.6, ar=c(-0.4), innov=innov.series, exact.innov=TRUE)
arfima.sim <- function(n, d=0, ar=numeric(), ma=numeric(), mu=0, sig2=1, stat.int=FALSE, 
                       n.burn, innov, exact.innov=TRUE){

  if (missing(n.burn)){
    if (d>=0) n.burn <- 100 else n.burn <- 10
  }
  if (missing(innov)) innov <- rnorm(n+n.burn,0,sig2)
    
  if (length(innov) < n + n.burn){
    # If given innovations aren't long enough for desired series length, either turn off burn-in or preprend innovations with samples from it
    # If given innovations are too long, this effectively lengthens the burn-in length. For d>0, this can increase variability in sample mean and variance of output.
    if (exact.innov){
      n.burn <- 0
    } else {
      n.diff <- n + n.burn - length(innov)
      innov <- c(sample(innov, n.diff, replace=TRUE), innov)
    }
  }
  
  if (stat.int && d>=0.5){
    m <- m0(d)
    d <- Phi0(d)
    # out <- better.filter(innov, arfima.to.ma(d, ar, ma, length(innov)-1))
    out <- fft.conv(diff.ma(d, length(innov)), innov)
    out <- fft.conv(arma.to.ma(ar, ma, length(out) - 1), out) + mu
    for (i in 1:m) out <- cumsum(out)
  } else {
    # out <- better.filter(innov, arfima.to.ma(d, ar, ma, length(innov)-1)) + mu
    m <- m0(d)
    d <- Phi0(d)
    out <- innov
    if (m>0){ 
      for (i in 1:m) out <- cumsum(out)
    } else if (m<0){
      out <- diff(out, differences = -m)
    }
    out <- fft.conv(diff.ma(d, length(out)), out)
    out <- fft.conv(arma.to.ma(ar, ma, length(out) - 1), out) + mu
  }
  
  # Take the most recent part of the sample
  out <- out[(length(out) - n + 1):length(out)]
  
  
  
  return(out)
}

