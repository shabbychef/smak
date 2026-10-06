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

# Created: 2026-09-22
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' Scalable Model Averaging.
#' 
#' @section Frequentist Model Averaging:
#'
#' We consider the basic problem of regressing a scalar \eqn{y} against a vector of
#' independent variables \eqn{\vec{x}}.
#' There is a long history of practitioners applying some form of \emph{variable selection} to select
#' a subset of the vector \eqn{\vec{x}} to use in the regression modeling,
#' as a way to avoid overfitting and improve generalization to out-of-sample data.
#' Besides the intractability of testing all \eqn{2^p} subsets of the \eqn{p} variables in \eqn{\vec{x}},
#' variable selection is often a precursor to questionable statistical practice.
#' 
#' As an alternative to variable selection, Hjort and Claeskens introduced 
#' \href{https://www.jstor.org/stable/30045339}{Frequentist Model Averaging}.
#' The idea behind FMA is to average together the regression coefficients
#' built on different subsets of the \eqn{\vec{x}}, with weights determined
#' typically by one of two methods, Mallows model averaging or Jackknife model averaging.
#'
#' @section Scalable Model Averaging:
#' 
#' While FMA supplies a stronger theoretical basis than variable selection,
#' it does not solve the computational issues.
#' Enter _Scalable Model Averaging_. 
#' Under this technique one first computes the singular value decomposition of the
#' matrix of observed independent variables, \eqn{X}.
#' Then one performs \eqn{p} regressions of the \eqn{y} against each of the left singular values of the \eqn{X},
#' then performs model averaging on only those \eqn{p} regressions.
#' Zhu \emph{et al.} show how to do this, and prove that the performance is equivalent
#' to the computationally intractable version run on all \eqn{2^p} regressions.
#'
#' This package implements that process of performing SVD then regression and
#' weight selection. The regressions produce model fits, but not estimates of
#' uncertainty of the coefficients or total model error. For that one would
#' have to rely on, for example, bootstrapping.
#'
#' @section Legal Mumbo Jumbo:
#'
#' smak is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-zhu
#'
#' @name smak
#' @rdname smak
#' @docType package
#' @title scalable model averagking kit
#' @keywords package
#' @note
#' 
#' This package is maintained as a hobby. 
#'
"_PACKAGE"

#' @title News for package 'smak':
#'
#' @description 
#'
#' News for package \sQuote{smak}
#'
#' \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\smak}{\CRANpkg{smak}}
#'
#' @section \smak{} Initial Version 0.1.0 (2026-10-06) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name smak-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:et:nu:syn=r:ft=r
