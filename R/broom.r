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
# Created: 2026.09.23
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @title broom .
#'
#' @description 
#'
#' One sentence or so that tells you some more.
#'
#' @details
#'
#' Really detailed. \eqn{\zeta}{zeta}.
#'
#' A list:
#' \itemize{
#' \item I use \eqn{n}{n} to stand for blah.
#' \item and so forth....
#' }
#' \describe{
#' \item{a}{value.}
#' \item{b}{factor.}
#' }
#'
#' @usage
#'
#' broom(x, n, zeta, ...)
#'
#' @param x vector of blah
#' @param n number of blah
#' @param ... arguments passed on to ...
#' @inheritParams dt
#' @inheritParams same_package_function
#' @inheritParams anotherPackage::function
#' @template param-ope
#' @return \code{dsr} gives the density, \code{psr} gives the distribution function,
#' \code{qsr} gives the quantile function, and \code{rsr} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @keywords distribution 
#' @keywords io
#' @keywords plotting
#' @aliases psr qsr rsr
#' @seealso t-distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' @note
#' This is a thin wrapper on the t distribution. 
#' @template etc
#' @template sr
#' @template smak
#' @references
#'
#' Johnson, N. L., and Welch, B. L. "Applications of the non-central t-distribution."
#' Biometrika 31, no. 3-4 (1940): 362-389. \url{http://dx.doi.org/10.1093/biomet/31.3-4.362}
#'
#' @examples 
#' y <- broom(20, 10)
#' \dontrun{
#' y <- broom(20, 10)
#' }
#' @author Steven E. Pav \email{steven@@gilgamath.com}
#' @export
broom <- function(x, n, zeta, type=c('mean','median','mad'),...) {
	type <- match.arg(type)



}

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics augment
#' @export
generics::augment

#' @title Tidy
#'
#' @description
#'
#' Tidy a smam model.
#'
#' @details
#'
#' Returns a table with information on the fit coefficients
#' from the scalable model averaging model.
#' At the moment, no standard errors are computed.
#'
#' @param x an object of type \code{smam}.
#' @param ... arguments for generic consistency.
#' @return A tidy \code{tibble::tibble()} with fields
#' \describe{
#'  \item{term}{The name of the estimated parameter. Betas typically come before gammas.}
#'  \item{estimate}{The estimated parameter.}
#' }
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
#' print(tidy(mod0))
#'
#' @template etc
#' @rdname tidy
#' @export
tidy.smam <- function(x, ...) {
	result <- data.frame(term=names(x$beta), estimate=as.numeric(x$beta))
	return(result)
}

#' @title Augment
#'
#' @description
#'
#' Augment a smam model.
#'
#' @details
#'
#' Returns a table with information on the overall fit
#' an estimated Harville or Henery model.
#'
#' @param x an object of type \code{harsm} or \code{hensm}
#' @param ... arguments for generic consistency.
#' @return A glanced \code{tibble::tibble()} with fields
#' \describe{
#'  \item{df}{The degrees of freedom of the model.}
#'  \item{logLik}{The log-likelihood of the model.}
#'  \item{AIC}{Akaike's Information Criterion for the model.}
#'  \item{nobs}{The number of observations, if this is available, otherwise ‘NA’.}
#' }
#' @note
#' In the future this may include information about the regularization, if any.
#' @seealso \code{\link[maxLik]{glance.maxLik}}.
#'
#' @examples
#'
#' # softmax on the Best Picture data
#' data(best_picture)
#' df <- best_picture
#' df$place <- ifelse(df$winner,1,2)
#' df$weight <- ifelse(df$winner,1,0)
#'
#' fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor + Drama
#' fit0 <- harsm(fmla,data=df,group=year,weights=weight)
#' print(glance(fit0))
#'
#' @template etc
#' @rdname glance
#' @export
augment.smam <- function(x, 
       data = model.frame(x),
       newdata = NULL,
       se_fit = FALSE,
			 ...) {
	stopifnot(!se_fit)
	# ack what to do here?
	#df <- broom:::augment_newdata(x, data, newdata, se_fit, ...)
	retv <- predict(x, newdata)
	result <- data.frame(retv)
	return(result)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
