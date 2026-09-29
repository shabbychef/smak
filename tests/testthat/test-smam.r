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

context("numerical robustness") #FOLDUP
test_that("mallows method works when a feature is orthogonal to response", {
  # avoid 1/0 error in inverse R matrix
  X <- diag(3)
  y <- c(1, 1, 0) # y is orthogonal to the 3rd singular vector
  expect_error(smamfit(y=y, X=X, method='mallows'), NA)
})
test_that("jackknife method handles saturated models / perfect leverage", {
  # avoid (1 - U^2) = 0 in jackknife
  X <- diag(3)
  y <- c(1, 2, 3)
  expect_error(smamfit(y=y, X=X, method='jackknife'), NA)
})
test_that("sirs throws appropriate error for n=2", {
  # avoid divide by 0 when n=2 by throwing error.
  X <- matrix(rnorm(10), nrow=2)
  y <- rnorm(2)
  expect_error(sirs(y, X), "n > 2") 
})
test_that("predict.smam handles factors with missing levels in newdata", {
  # predict even when levels are missing
  df <- data.frame(y=rnorm(4), f=factor(c("A", "B", "A", "B")), x=rnorm(4))
  mod <- smam(y ~ f + x, data=df)
  df_new <- data.frame(f=factor(c("A")), x=rnorm(1)) # "B" is missing
  expect_error(predict(mod, newdata=df_new), NA)
})
test_that("smamfit correctly handles small alpha values", {
  # when alpha is tiny, do not return only zero columns. 
  X <- matrix(rnorm(30), nrow=3)
  y <- rnorm(3)
  expect_error(smamfit(y=y, X=X, method='mallows', 
                       control=list(wide_pragma="alpha", alpha=1e-10)), NA)
})
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
test_that("model and predict handle offset", { #FOLDUP
  nobs <- 100
  nfeat <- 5
  set.seed(123)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
  # Add some known offset
  off_set <- rnorm(length(eta))
  y <- rnorm(length(eta), mean = eta + off_set)

  Xdf <- as.data.frame(X)
  varnames <- colnames(Xdf)
  Xdf$y <- y
  Xdf$off_set <- off_set

  fmla <- as.formula(paste0("y ~ -1 + offset(off_set) + ", paste(varnames, collapse=" + ")))
  expect_error(afit <- smam(fmla, Xdf), NA)

  # fitted values should be roughly y, and they should match predict()
  prd1 <- predict(afit)
  expect_equal(prd1, afit$fitted.values)

  # predict on new data with new offset
  newX <- matrix(rnorm(10 * nfeat), ncol = nfeat)
  newXdf <- as.data.frame(newX)
  colnames(newXdf) <- varnames
  newXdf$off_set <- rnorm(10)
  expect_error(prd2 <- predict(afit, newdata = newXdf), NA)
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
#UNFOLD

context("wide data") #FOLDUP
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
  expect_error(afit <- smamfit(y, X, control=list(wide_pragma='SIRS')),NA)
  expect_error(afit <- smamfit(y, X, control=list(wide_pragma='SIRSu')),NA)
  # throw.
  expect_error(afit <- smamfit(y, X, control=list(wide_pragma='error')))
  expect_error(afit <- smamfit(y, X, control=list(wide_pragma='unknown method')))
}) #UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
