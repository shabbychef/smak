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
# Created: 2026.09.25
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# Computes the \hat{w} vector to perform SIRS variable selection.
.what <- function(y, X, ties_pragma=c('ignore','half_value')) {
  ties_pragma <- match.arg(ties_pragma)
  sidx <- sort(y, index.return=TRUE)
  y <- sidx$x
  X <- X[sidx$ix, ,drop=FALSE]
  n <- nrow(X)
  switch(ties_pragma,
  ignore={
    wtilde <- rep(0, ncol(X))
    for (jidx in (2:n)) {
      iidx <- which(y[1:(jidx-1)] < y[jidx])
      wtilde <- wtilde + (colSums(X[iidx,,drop=FALSE]))^2
    }
  },
  half_value={
    wtilde <- rep(0, ncol(X))
    for (jidx in (2:n)) {
      less_idx <- which(y[1:(jidx-1)] < y[jidx])
      eq_idx <- which(y == y[jidx])
      wtilde <- wtilde + (colSums(X[less_idx,,drop=FALSE]) + 0.5*(colSums(X[eq_idx,,drop=FALSE]) - X[jidx,]))^2
    }
  })
  # this is a U statistic
  what <- wtilde / (n*(n-1)*(n-2))
  return(what)
}

# Performs the SIRS algorithm of zhu et al.
# Modified to deal with ties in the y.
# @template ref-sirs

#' @title Model-Free Feature Screening
#'
#' @description 
#'
#' Performs the model-free feature screening procedure of Zhu et al. for
#' regression problems with many features.
#'
#' @details
#'
#' Performs a model-free screening procedure. 
#' First each column of \eqn{X} is z-scored to get zero mean unit variance
#' variables. Then the \eqn{\hat{w}} is computed for each of the columns of
#' \eqn{X} against \eqn{y}. The kth element of this vector is computed as
#' \deqn{\frac{1}{n(n-1)(n-2)}sum_{1\le j \le n}\left(\sum_{1\le i \le n}
#' X_{i,k} I(y_i < y_j)\right)^2.}
#' Then we select indices \eqn{k} based on the values of this \eqn{\hat{w}},
#' based on the thresholding:
#'
#' \describe{
#' \item{hard}{Select a fixed number of variates, based on the ordering of
#' \eqn{\hat{w}}.}
#' \item{soft}{Variates are selected if they are larger than the maximum of
#' \eqn{\hat{w}} computed on purely noise variates.}
#' }
#'
#' @param thresholding  the method for performing thresholding of the
#' \eqn{\hat{w}} variable to select predictors. For soft thresholding we
#' oversample with noise; for hard thresholding we pick the top \eqn{N} values.
#' @param N  the number of variables to select for hard thresholding. Defaults
#' to \eqn{n / log(n)}.
#' @param dprop  for soft thresholding controls the number of noise variables
#' added to the set. We add \eqn{dprop p} variables.
#' @param ties_pragma  the method for dealing with ties in the \eqn{y} values.
#' For \sQuote{ignore} we ignore them and implement the classic value of
#' Zhu et. al.  For \sQuote{half_value} the indicator function returns one half
#' for tied values in the \eqn{y}.
#' @inheritParams smamfit 
#' @return A sorted vector of the indices of which columns of \code{X} would be
#' selected, with values in \eqn{1, 2, \ldots, p}.
#' @seealso \code{\link{smamfit}}.
#' @template etc
#' @template ref-sirs
#' @examples 
#' set.seed(1234)
#' X <- matrix(rnorm(200*500),ncol=500)
#' beta <- 2 * round(rnorm(ncol(X), sd=0.2)) * runif(ncol(X))
#' sum(abs(beta) > 0)
#' true_vars <- which(abs(beta) > 0)
#' y <- rnorm(nrow(X), mean=X %*% beta)
#' soft_vals <- sirs(y, X, thresholding='soft')
#' all(soft_vals %in% true_vars)
#' 
#' hard_vals <- sirs(y, X, thresholding='hard')
#' all(true_vars %in% hard_vals)
#' @export
sirs <- function(y, X, thresholding=c('hard','soft'), N=NULL, d_prop=1.0,
                 ties_pragma=c('ignore', 'half_value')) {
  thresholding <- match.arg(thresholding)
  ties_pragma <- match.arg(ties_pragma)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(n == length(y))
  stopifnot(n > 2)
  
  # Standardize X to have zero mean and unit variance
  colsd <- apply(X, 2, FUN=sd)
  colsd[colsd == 0] <- 1
  X <- scale(X, center=TRUE, scale=colsd)

  switch(thresholding,
         hard={ 
           if (is.null(N)) {
             N <- ceiling(n/log(n))
           }
           what <- .what(y, X, ties_pragma=ties_pragma)
           sort_what <- sort(what, decreasing=TRUE, index.return=TRUE)
           retval <- sort_what$ix[1:N]
         },
         soft={ 
           n_noise <- ceiling(d_prop * p)
           what <- .what(y, cbind(X, matrix(rnorm(n_noise*n),nrow=n)), ties_pragma=ties_pragma)
           C_d <- max(what[p+(1:n_noise)])
           retval <- which(what[1:p] > C_d)
         }
  )
  retval <- sort(retval)
  return(retval)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
