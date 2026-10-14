dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Scalable Model Averageing Kit
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Supports frequentist model averaging on linear
    regression problems via the scalable method of
    Zhu et. al. (2022) <doi:10.1080/07350015.2022.2116442>.
    The scalable model averaging technique projects to
    singular vectors, performs one linear regression per singular
    vector, then combines the estimates to optimize either
    Mallows' criterion, or the jackknife estimate of error.
Depends: 
    R (>= 4.0)
Imports:
    generics,
    methods,
    quadprog
dnl LinkingTo: Rcpp
Suggests: 
    broom,
    testthat, 
    knitr
URL: https://github.com/shabbychef/PKG_NAME()
dnl VignetteBuilder: knitr
Encoding: UTF-8
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4:et
