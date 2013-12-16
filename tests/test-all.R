if (require(testthat)) {
  library(lava)
  if (exists("test_check"))
      test_check("lava")
}

