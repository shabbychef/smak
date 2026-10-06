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
test_that("coef", { #FOLDUP
  nfeat <- 5
  nobs <- 100
  set.seed(1234)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
  y <- rnorm(length(eta), mean=eta)
  expect_error(afit <- smamfit(y, X), NA)
  expect_error(resp <- coef(afit), NA)
  expect_equal(as.numeric(afit$beta), as.numeric(resp))
  expect_error(resp <- coefficients(afit), NA)
  expect_equal(as.numeric(afit$beta), as.numeric(resp))
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

context("vs OLS") #FOLDUP
test_that("sma gives results like OLS for large n", {
  nobs <- 100000
  nprd <- 3
  set.seed(1234)
  X <- matrix(rnorm(nobs*nprd), ncol=nprd)
  beta <- runif(nprd,min=1,max=2)
  y <- rnorm(nobs,mean=X %*% beta)
  expect_error(mmafit <- smamfit(y=y, X=X, method='mallows'),NA)
  expect_error(jmafit <- smamfit(y=y, X=X, method='jackknife'),NA)
  expect_error(olsfit <- lm.fit(x=X, y=y), NA)
  expect_equal(as.numeric(mmafit$beta), as.numeric(olsfit$coefficients), tolerance=1e-4)
  expect_equal(as.numeric(jmafit$beta), as.numeric(olsfit$coefficients), tolerance=1e-4)
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
test_that("fitted.values OK", { #FOLDUP
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
  expect_error(bfit <- smam(fmla, Xdf), NA)
  expect_error(prd2 <- predict(bfit, newdata=Xdf), NA)
	expect_equal(prd2, bfit$fitted.values)

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
test_that("model and predict handle factors", { #FOLDUP
  nobs <- 100
  nfeat <- 5
	nlevl <- 3
  set.seed(456)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  eta <- X %*% beta
	# now do something here with the factor ...
	ebeta <- runif(nlevl)
	fac_idx <- sample(seq_along(ebeta), size=nobs, replace=TRUE)
	eta <- eta + ebeta[fac_idx]
  y <- rnorm(length(eta), mean = eta)

  Xdf <- as.data.frame(X)
  varnames <- colnames(Xdf)
  Xdf$y <- y
  Xdf$fac <- factor(letters[fac_idx])

  fmla <- as.formula(paste0("y ~ -1 + fac + ", paste(varnames, collapse=" + ")))
  expect_error(afit <- smam(fmla, Xdf), NA)

  # fitted values should be roughly y, and they should match predict()
  prd1 <- predict(afit)
  expect_equal(prd1, afit$fitted.values)

  # predict on new data with new offset
  newX <- matrix(rnorm(10 * nfeat), ncol = nfeat)
  newXdf <- as.data.frame(newX)
  colnames(newXdf) <- varnames
	for (idx in seq_len(nlevl)) {
		newXdf$fac <- factor(letters[idx], levels=letters[seq_len(nlevl)])
		expect_error(prd2 <- predict(afit, newdata = newXdf), NA)
	}
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

context("weights") #FOLDUP
test_that("smamfit weights act as replication weights", { #FOLDUP
  nfeat <- 4
  nobs <- 40
  set.seed(9876)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  y <- rnorm(nobs, mean = X %*% beta)
  wt <- rep(2, nobs)

  idx <- rep(seq_len(nobs), times = wt)
  X_dup <- X[idx, ]
  y_dup <- y[idx]

  for (method in c('jackknife', 'mallows')) {
    expect_error(afit <- smamfit(y, X, wt = wt, method = method), NA)
    expect_error(bfit <- smamfit(y_dup, X_dup, method = method), NA)
    expect_equal(as.numeric(afit$beta), as.numeric(bfit$beta), tolerance = 1e-6)
    expect_equal(as.numeric(afit$weights), as.numeric(bfit$weights), tolerance = 1e-6)
    expect_equal(as.numeric(afit$sigma2), as.numeric(bfit$sigma2), tolerance = 1e-6)
  }

  wt <- sample(1:4, nobs, replace = TRUE)

  idx <- rep(seq_len(nobs), times = wt)
  X_dup <- X[idx, ]
  y_dup <- y[idx]

  for (method in c('jackknife', 'mallows')) {
    expect_error(afit <- smamfit(y, X, wt = wt, method = method), NA)
    expect_error(bfit <- smamfit(y_dup, X_dup, method = method), NA)
    expect_equal(as.numeric(afit$beta), as.numeric(bfit$beta), tolerance = 1e-6)
    expect_equal(as.numeric(afit$weights), as.numeric(bfit$weights), tolerance = 1e-6)
    expect_equal(as.numeric(afit$sigma2), as.numeric(bfit$sigma2), tolerance = 1e-6)
  }

}) #UNFOLD
test_that("smam formula weights act as replication weights", { #FOLDUP
  nfeat <- 4
  nobs <- 40
  set.seed(5432)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  y <- rnorm(nobs, mean = X %*% beta)
  wts <- sample(1:4, nobs, replace = TRUE)

  df <- as.data.frame(X)
  df$y <- y
  df$wts <- wts
  fmla <- as.formula(paste0("y ~ ", paste(names(as.data.frame(X)), collapse = " + ")))

  idx <- rep(seq_len(nobs), times = wts)
  df_dup <- df[idx, ]
  df_dup$wts <- NULL

  # Test bare symbol, string, and vector weight specifications
  expect_error(fit_sym <- smam(fmla, df, weights = wts), NA)
  expect_error(fit_str <- smam(fmla, df, weights = "wts"), NA)
  expect_error(fit_vec <- smam(fmla, df, weights = df$wts), NA)
  expect_error(fit_dup <- smam(fmla, df_dup), NA)

  expect_equal(as.numeric(fit_sym$beta), as.numeric(fit_str$beta))
  expect_equal(as.numeric(fit_sym$beta), as.numeric(fit_vec$beta))
  expect_equal(as.numeric(fit_sym$beta), as.numeric(fit_dup$beta), tolerance = 1e-6)
  expect_equal(as.numeric(fit_sym$weights), as.numeric(fit_dup$weights), tolerance = 1e-6)
}) #UNFOLD
test_that("zero weights drop observations like row deletion", { #FOLDUP
  nfeat <- 4
  nobs <- 40
  set.seed(1357)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  y <- rnorm(nobs, mean = X %*% beta)
  wt <- sample(1:4, nobs, replace = TRUE)
  wt[c(3, 7, 12)] <- 0

  idx <- rep(seq_len(nobs), times = wt)
  X_dup <- X[idx, ]
  y_dup <- y[idx]

  for (method in c('jackknife', 'mallows')) {
    expect_error(afit <- smamfit(y, X, wt = wt, method = method), NA)
    expect_error(bfit <- smamfit(y_dup, X_dup, method = method), NA)
    expect_equal(as.numeric(afit$beta), as.numeric(bfit$beta), tolerance = 1e-6)
  }
}) #UNFOLD
test_that("high weights converge to OLS", { #FOLDUP
  nfeat <- 3
  nobs <- 100
  set.seed(1212)
  X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
  beta <- rnorm(nfeat)
  y <- rnorm(nobs, mean = X %*% beta)
  wt <- rep(10000,nobs)

  expect_error(olsfit <- lm.fit(x=X, y=y), NA)
  for (method in c('jackknife', 'mallows')) {
    expect_error(bfit <- smamfit(y, X, wt = wt, method = method), NA)
    expect_equal(as.numeric(olsfit$coefficients), as.numeric(bfit$beta), tolerance = 1e-4)
    #expect_error(afit <- smamfit(y, X, method = method), NA)
    #expect_equal(as.numeric(afit$beta), as.numeric(bfit$beta), tolerance = 1e-6)
  }
}) #UNFOLD
test_that("very large weights do not crash the fit", { #FOLDUP
  nfeat <- 5
  nobs <- 100
  set.seed(345)
	replicate(5,{
    X <- matrix(rnorm(nobs * nfeat), ncol = nfeat)
    y <- rnorm(nobs)
    wt_large <- 10^runif(nobs, 3, 9)
    for (method in c('jackknife', 'mallows')) {
      expect_error(afit <- smamfit(y, X, method=method, wt = wt_large), NA)
    }
  })
}) #UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
