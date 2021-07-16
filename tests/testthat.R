Sys.setenv("R_TESTS" = "")

library(testthat)
library(scMethrix)

test_check("scMethrix")
