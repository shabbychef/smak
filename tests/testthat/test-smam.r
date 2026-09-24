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

test_that("smamfit", {
  #FOLDUP
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
}) #UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
