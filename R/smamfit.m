# /usr/bin/r
#
# Copyright 2026-2026 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of smak.
#
# smak is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# smak is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with smak.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2026.09.22
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# apparently necessary to register smam as an S3 class
setOldClass('smam')

#' @title Auxiliary for Controlling SMA Fitting
#'
#' @description 
#'
#' Auxiliary function for \code{smam} fitting.
#'
#' @details
#'
#' The \code{control} argument of \code{smam} is by default
#' passed to the \code{control} argument of \code{smamfit}, which uses its
#' selements as arguments to \code{smam.control}: the latter provides
#' defaults and sanity checking.
#'
#' The wide pragma controls how fitting deals with wide matrices \eqn{X}. 
#' \describe{
#'  \item{SIRS}{Implements the SIRS procedure of Zhu et al. on \eqn{X}.}
#'  \item{SIRSu}{Implements the SIRS procedure on the left singular values
#'  \eqn{U}.}
#'  \item{alpha}{Performs SVD, then keeps the singular values that account for
#'  at least \eqn{\alpha} of the total sum of singular values.}
#'  \item{fixed_k}{Keep a fixed number of singular values.}
#'  \item{error}{Throws an error.}
#' }
#'
#' The tall pragma controls how fitting deals with tall matrices \eqn{X}.
#' \describe{
#'  \item{none}{Takes all singular values.}
#'  \item{tolerance}{Keeps singular values up to machine epsilon ratio to the
#'  largest singular value.}
#' }
#'
#' @param wide_pragma  what to do if X has more columns than rows, see Details.
#' @param tall_pragma  what to do if X has as many rows as columns or more, see
#' Details.
#' @param alpha  the \eqn{alpha} for the \dQuote{alpha} wide pragma.
#' @param k  the number of singular values kept for the \dQuote{fixed_k} wide
#' pragma.
#' @return A list.
#' @template ref-zhu
#' @template ref-sirs
#' @seealso \code{\link{smamfit}}
#' @template etc
#' @export
smam.control <- function(wide_pragma=c("SIRS", "SIRSu", "alpha", "fixed_k", "error"),
												tall_pragma=c("tolerance", "none"),
												alpha=0.99,
												k=100) {
	wide_pragma <- match.arg(wide_pragma)
	list(wide_pragma=wide_pragma, alpha=alpha, k=k)
}

.do_svd <- function(X, control=list(), ...) {
	control <- do.call("smam.control", control)
	if (nrow(X) < ncol(X)) {
		switch(control$wide_pragma,
					 SIRS={ stop("not yet implemented") },
					 SIRSu={ 
						 USV <- svd(X)
						 stop("not yet implemented") 
					 },
					 fixed_k={
						 nuv <- min(nrow(X), control$k)
						 USV <- svd(X, nu=nuv, nv=nuv)
					 },
					 alpha={
						 USV <- svd(X)
						 k <- sum(cumsum(USV$d) < control$alpha * sum(USV$d))
						 USV$u <- USV$u[, 1:k, drop = FALSE]
						 USV$d <- USV$d[1:k]
						 USV$v <- USV$v[, 1:k, drop = FALSE]
					 },
					 error={ stop("Will not perform SVD on wide matrix under this pragma.") })
	} else {
		# perform SVD
		USV <- svd(X)
		switch(control$tall_pragma,
			none={
				# noop
			},
			tolerance={
				tol <- max(n, p) * USV$d[1] * .Machine$double.eps
				k <- sum(USV$d > tol)
				USV$u <- USV$u[, 1:k, drop = FALSE]
				USV$d <- USV$d[1:k]
				USV$v <- USV$v[, 1:k, drop = FALSE]
			}
		)
	}
	return(USV)
}



#' @title scalable frequentist model averaging regression fit.
#'
#' @description 
#'
#' Fits the scalable model averaging model to a regression problem.
#'
#' @details
#'
#' Finds a vector of coefficients \eqn{\beta} such that \eqn{X\beta} is 
#' approximately equal to \eqn{y}. Does this by computing a singular
#' value decomposition of \eqn{X}, then performs \eqn{p} linear regressions
#' of the column of the left singular vectors against \eqn{y}, then
#' takes a weighted combination of those. 
#'
#' Supports both the Mallows Model Averaging and Jackknife Model Averaging
#' approaches.
#'
#' @param y vector of response data, of length \eqn{n}.
#' @param X a \eqn{n \times p}{n x p} matrix of independent variables.
#' @param method  the selection method for the weights, either Jackknife model
#' averaging or Mallows' model averaging.
#' @param wt an optional vector of length \eqn{n} giving the weight of each
#' observation. If given must be non-negative. Weights are replication weights.
#' @param sigma2 an optional squared variance to be used in Mallows' model
#' averaging. If not given, estimated from the data.
#' @param Xnames an optional vector of the names of the columns of X. If null
#' or not given, defaults to \code{colnames(X)}. Used only for naming the
#' output.
#' @param return.U  a Boolean for whether the matrix U should be returned in
#' the output. 
#' @param return.X  a Boolean for whether the matrix X should be returned in
#' the output. 
#' @param return.fitted  a Boolean for whether the yhat fitted values should be returned in
#' the output. 
#' @param control a list of parameters for controlling the the fitting process.
#' This is passed to \code{\link{smam.control}}.
#' @return An object of class \code{smam}.
#' @keywords fitting
#' @seealso the friendly interface \code{\link{smam}}.
#' @template etc
#' @template smak
#' @template ref-zhu
#' @examples 
#'
#' nobs <- 100
#' nfeat <- 5
#' set.seed(1234)
#' X <- matrix(rnorm(nobs * nfeat),ncol=nfeat)
#' beta <- rnorm(nfeat)
#' eta <- X %*% beta
#' y <- rnorm(ncol(X), mean=X %*% beta, sd=1)
#'
#' mod0 <- smamfit(y=y,X=X)
#' summary(mod0)
#' @export
smamfit <- function(y, X, method=c('jackknife','mallows'), wt=NULL, sigma2=NULL, 
										Xnames=NULL, return.U=FALSE, return.X=TRUE, return.fitted=TRUE, 
										control=list(...), ...) {
	method <- match.arg(method)

  y <- as.numeric(y)
	n <- nrow(X)
	p <- ncol(X)
	stopifnot(n == length(y))

	if (is.null(Xnames)) {
		Xnames <- colnames(X)
		if (is.null(Xnames)) {
			Xnames <- paste0("X_",seq_len(p))
		}
	}

	if (!is.null(wt)) {
		stopifnot(n == length(wt), all(wt >= 0))
		root_w <- sqrt(wt)
    y <- root_w * y
		USV <- .do_svd(root_w * X, control=control)
	} else {
		USV <- .do_svd(X, control=control)
	}


	U <- USV$u
  D_inv <- 1 / USV$d
	V <- USV$v

	k <- ncol(U)
	Unames <- colnames(U)
	if (is.null(Unames)) {
		Unames <- paste0("U_",seq_len(k))
	}

	# compute the univariate betas, which is are the betas of each U against y
	betas <- t(U) %*% y
	# compute the mu hats
	muhats = t(t(U) * as.numeric(betas))
	# compute the weights hat{w}_j
	switch(method,
		jackknife={
			Dmat <- 2 * crossprod(t(t(muhats) - y) / (1 - U^2))
			dvec <- rep(0, k)
			factorized <- FALSE
			# Ensure Dmat is positive definite for solve.QP numerical stability
			diag(Dmat) <- diag(Dmat) + 1e-10
		},
		mallows={
			# dMatrix will be 2 * diag(betas^2), and we want the inverse of the
			# square root of this, so (1/sqrt(2)) * diag(betas^-1)
			invRmat <- (1/sqrt(2)) * diag(1 / as.numeric(betas), nrow=k)
			b2 <- betas^2
			if (is.null(sigma2)) {
				rss_full <- sum(y^2) - sum(b2)
				if (n > k) {
					sigma2 <- rss_full / (n - k)
				} else {
					# Fallback for saturated cases
					sigma2 <- rss_full / n
				}
			} 
			Dmat <- invRmat
			dvec <- 2 * (b2 - sigma2)
			factorized <- TRUE
		})

	# constraints are sum(w) = 1 and w_j >= 0
  Amat <- cbind(rep(1, k), diag(k))
  bvec <- c(1, rep(0, k))
	meq <- 1
	
	# solve the quadratic program
  res   <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq, factorized=factorized)

	# interpret
  w_opt <- res$solution
  
  # Numerical cleanup
  w_opt[w_opt < 0] <- 0
  w_opt <- w_opt / sum(w_opt)
  
  # Fitted values and beta in original (unscaled) space.
  # In scaled space: fitted_tilde = U %*% (w_opt * a) = W^{1/2} fitted_y.
  # Back-transform by dividing by root_w when observation weights were supplied.
  theta_hat <- betas * w_opt
  
  # beta = V D^{-1} theta_hat  (gives WLS estimate in original parameterization)
  beta_hat <- as.vector(V %*% (D_inv * theta_hat))
  
  retv <- list(
    beta = beta_hat,
    weights = w_opt,
    univariate_betas = betas,
    sigma2 = sigma2
  )
	names(retv$beta) <- Xnames
	names(retv$weights) <- Unames
	names(retv$univariate_betas) <- Unames
	if (return.U) { retv$U <- U }
	if (return.X) { retv$X <- X }
	if (return.fitted) {
		fitted_scaled <- as.vector(U %*% theta_hat)
		retv$fitted.values <- if (!is.null(wt)) fitted_scaled / root_w else fitted_scaled
	}
	class(retv) <- 'smam'
	return(retv)
}

#' @title Friendly interface to scalable model averaging regression.
#'
#' @description
#'
#' A user friendly interface to the scalable frequentist model averaging
#' regression procedure.
#'
#' @details
#'
#' Performs the scalable model averaging procedure to estimate linear
#' regression of \eqn{y} against \eqn{X}.
#' We find vector of coefficients \eqn{\beta} such that \eqn{X\beta} is 
#' approximately equal to \eqn{y}.
#'
#' @inheritParams stats::lm
#' @inheritParams smamfit   
#' @param na.action  How to deal with missing values in the outcomes,
#' weights, etc.
#' @param weights  an optional vector of weights, or the string or bare name of the
#' weights in the \code{data} for use in the fitting process. Set to \code{NULL}
#' for none.
#' @template etc
#' @return An object of class \code{smam}.
#' @keywords fitting
#' @seealso \code{\link{smamfit}}.
#'
#' @examples
#'
#' nobs <- 100
#' nfeat <- 5
#' set.seed(1234)
#' X <- matrix(rnorm(nobs * nfeat),ncol=nfeat)
#' beta <- rnorm(nfeat)
#' y <- rnorm(ncol(X), mean=X %*% beta, sd=1)
#' # now the pretty frontend
#' data <- cbind(data.frame(outcome=y),as.data.frame(X))
#'
#' fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
#' fitm <- smam(fmla,data)
#'
#' # with weights
#' data$wts <- runif(nrow(data),min=1,max=2)
#' fitm <- smam(fmla,data,group=race,weights=wts)
#'
#' @importFrom stats coef formula model.frame model.matrix na.omit model.response model.weights
#' @export
#' @rdname smam
smam <- function(formula,data,weights=NULL,na.action=na.omit,method=c('jackknife','mallows'), ...) {
	substitute(formula)
	# I find it highly offensive that this cannot be done reasonably
	# easily in a subfunction because of NSE whatever.

	# https://stackoverflow.com/q/53827563/164611
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	strmf <- as.character(mf)
	#turn weights into symbol if character is passed
	if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
	m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE 
	mf[[1L]] <- quote(stats::model.frame) 
	mf <- eval(mf, parent.frame()) #evaluate call

	X <- model.matrix(formula,mf)
	y <- as.vector(model.response(mf))
	wt <- as.vector(model.weights(mf))

	dat <- list(X=X,y=y,wt=wt,Xnames=colnames(X))
	# call the fit function
	retv <- smamfit(y=dat$y, X=dat$X, wt=dat$wt, method=method, Xnames=dat$Xnames, ...)
	retv$call <- cl
	retv$formula <- formula
	return(retv)
}


#' @rdname smam    
#' @importFrom stats predict
#' @param newdata  a \code{data.frame} from which we can extract a model
#' frame via the formula of the \code{object}.
#' @template param-group
#' @param ... other arguments.
#' @param type  indicates which prediction should be returned:
#' \describe{
#' \item{\code{eta}}{The odds.}
#' \item{\code{mu}}{The probability.}
#' \item{\code{erank}}{The expected rank.}
#' }
#' @param na.action  How to deal with missing values in \code{y}, \code{g},
#' \code{X}, \code{wt}, \code{eta0}.
#' @seealso \code{\link{smax}}, \code{\link{harsm_invlink}}.
#' @importFrom stats delete.response terms model.offset model.matrix model.extract as.formula na.pass
#' @export
#' @method predict smam
predict.smam <- function(
  object,
  newdata,
  na.action = na.pass,
  ...
) {
  fmla <- object$formula
  tt <- terms(fmla)
  Terms <- delete.response(tt)

  # https://stackoverflow.com/q/53827563/164611
  mf <- match.call(expand.dots = FALSE)
  # wheee
  names(mf) <- gsub('^newdata$', 'data', names(mf))

  m <- match(c("data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, 1L, m)]
  mf$drop.unused.levels <- TRUE
  # need this
  # mf$xlev <-
  mf[[1L]] <- quote(stats::model.frame)
  mf[[2L]] <- Terms
  mf <- eval(mf, parent.frame()) #evaluate call
  X <- model.matrix(as.formula(Terms), mf)
	# ooff, what do we do with this?
  y0 <- model.offset(mf)
  wt <- as.vector(model.weights(mf))

  dat <- list(X = X, y0 = y0, wt = wt)
  if (!all(colnames(dat$X) %in% names(object$beta))) {
    stop("some levels in data unknown to fit model")
  }
  # 2FIX: we subset beta to the colnames, but are these promised to be in the right order?
  yhat <- as.numeric(dat$X %*% matrix(object$beta[colnames(dat$X)], ncol = 1))
  # deal with offset
  if (!is.null(y0)) {
    yhat <- yhat + y0
  }
  attr(yhat, 'na.action') <- attr(dat, 'na.action')
  return(yhat)
}

#' @export
#' @importFrom stats printCoefmat
#' @importFrom methods show
#' @inheritParams base::print
#' @rdname smam
#' @method print smam
print.smam <- function(x, ...) {
	cat("Call:",'\n')
	show(x$call)
	cat('\n')
	cat("Coefficients:",'\n')
	show(x$beta)
	invisible(x)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
