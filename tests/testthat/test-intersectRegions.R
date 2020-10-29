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



# interactions ------------------------------------------------------------

context("bins")
test_that("bins1",{

  test_bw = GRanges(seqnames = "chr1",
                   ranges = IRanges::IRanges(
                     start = seq(1,100,by = 10),
                     width = 10
                   ),
                   score = 10:1)


  binSignal_cumSum(gr = test_bw,
                   min_value = 2,
                   n_bins = 2) -> result_bins

  expect_equal(max(result_bins$cumWidth,na.rm = T), 80)
  expect_equal(min(result_bins$cumWidth,na.rm = T), 10)
  expect_equal(sum(is.na(result_bins$cumWidth)),2)

  total_s1 = sum(width(result_bins[result_bins$bin_values == "eqFreqBin2of2"]))
  total_s2 = sum(width(result_bins[result_bins$bin_values == "eqFreqBin1of2"]))

  expect_equal(total_s1,total_s2)

})
