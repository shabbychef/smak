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
Description: utilities for performing regression tasks with 
    Frequentist model averaging, via the scalable method of
    Zhu et. al. (2022) <doi:10.1080/07350015.2022.2116442>.
Depends: 
    R (>= 4.0)
Imports:
    quadprog
dnl LinkingTo: Rcpp
Suggests: 
    testthat, 
    knitr
URL: https://github.com/shabbychef/PKG_NAME()
VignetteBuilder: knitr
Encoding: UTF-8
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4:et
