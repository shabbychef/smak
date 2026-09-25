# Copyright 2026 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

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

# env var:
# nb:
# see also:
# todo:
# changelog:
#
# Created: 2026.09.24
# Copyright: Steven E. Pav, 2026-2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
  set.seed(as.integer(charToRaw(str)))
}
#UNFOLD

context("smoke test") #FOLDUP
test_that("smamfit and smam", { #FOLDUP
  # travis only?
  #skip_on_cran()
  nfeat <- 5
	nobs <- 100
  set.seed(1234)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	expect_error(afit <- smamfit(y, X), NA)
	# same results from data frame
	Xdf <- as.data.frame(X)
	varnames <- names(afit$beta)
	colnames(Xdf) <- varnames
	Xdf$y <- y
	fmla <- as.formula(paste0("y ~ -1 + ",paste(varnames,collapse=" + ")))
	expect_error(bfit <- smam(fmla, Xdf), NA)
	expect_equal(afit$beta, bfit$beta)
}) #UNFOLD
test_that("mma and jma", { #FOLDUP
  nfeat <- 5
	nobs <- 100
  set.seed(4567)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	for (method in c('jackknife','mallows')) {
		expect_error(afit <- smamfit(y, X, method=method), NA)
	}
}) #UNFOLD
test_that("broom::tidy", { #FOLDUP
  nfeat <- 5
	nobs <- 100
  set.seed(1234)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	expect_error(afit <- smamfit(y, X), NA)
	expect_error(resp <- tidy(afit), NA)
	expect_equal(as.numeric(afit$beta), resp$estimate)
}) #UNFOLD
#UNFOLD

context("predict") #FOLDUP
test_that("predict method", { #FOLDUP
  nfeat <- 5
	nobs <- 100
  set.seed(4567)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	expect_error(afit <- smamfit(y, X), NA)
	expect_error(prd1 <- predict(afit), NA)
	expect_equal(prd1, afit$fitted.values)
	# now on new data.
  newX <- matrix(rnorm(10 * nfeat), ncol = nfeat)
	varnames <- names(afit$beta)
	newXdf <- as.data.frame(newX)
	colnames(newXdf) <- varnames
	# fails without a formula.
	expect_error(prd1 <- predict(afit, newdata=newXdf))

	# same results from data frame
	Xdf <- as.data.frame(X)
	varnames <- names(afit$beta)
	colnames(Xdf) <- varnames
	Xdf$y <- y
	fmla <- as.formula(paste0("y ~ -1 + ",paste(varnames,collapse=" + ")))
	expect_error(bfit <- smam(fmla, Xdf), NA)
	expect_error(prd2 <- predict(bfit, newdata=newXdf), NA)

	expect_error(afit2 <- smamfit(y, X, formula=fmla), NA)
	expect_error(prd3 <- predict(afit2, newdata=newXdf), NA)

}) #UNFOLD
test_that("broom::augment", { #FOLDUP
  nfeat <- 5
	nobs <- 100
  set.seed(4567)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	Xdf <- as.data.frame(X)
	varnames <- colnames(Xdf)
	Xdf$y <- y
	fmla <- as.formula(paste0("y ~ -1 + ",paste(varnames,collapse=" + ")))

	expect_error(afit <- smam(fmla, Xdf), NA)
	expect_error(prd1 <- augment(afit, newdata=Xdf), NA)
	# this fails.
	expect_error(prd2 <- augment(afit), NA)
}) #UNFOLD

context("wide data") 
test_that("wide methods", { #FOLDUP
	nobs <- 100
  nfeat <- nobs + 20
  set.seed(4567)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	y <- rnorm(length(eta), mean=eta)
	# fine
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='alpha',alpha=0.5)),NA)
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='fixed_k',k=10)),NA)
	# throw.
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='error')))
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='unknown method')))
	# not yet implemented.
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='SIRS')),'not yet implemented')
	expect_error(afit <- smamfit(y, X, control=list(wide_pragma='SIRSu')),'not yet implemented')
}) #UNFOLD


#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
