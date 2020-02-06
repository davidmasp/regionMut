test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})





# UTILS -------------------------------------------------------------------
context("utils")
test_that("utils1",{
  list("a" = c(2,3,4),
       "b" = c(1,2,2),
       "c" = c(2)) -> test

  regionmut_unlist(x = test) -> res

  expect_true(length(res) == length(unlist(test)))
  expect_equal(unique(names(res)), names(test))
})

