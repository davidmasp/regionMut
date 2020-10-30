# Bining


binSignal_uniqueRes <- function(gr,
                                    min_instances = 10,
                                    min_value = 1,
                                    n_bins = 4) {

  ## here I obtain the available resolutions of the original file.
  params_resolutions = table(GenomicRanges::width(gr))

  # this tells us how minimun number of regions are we allow to drop if
  # they are not at the same res as the rest of the fail. If it is bigger
  # than this number the function will break
  params_minInstances = min_instances

  if (any(params_resolutions < params_minInstances)){
    params_res <- names(params_resolutions[params_resolutions > params_minInstances])
  } else if (length(params_resolutions) == 1){
    params_res <- names(params_resolutions)
  } else {
    stop("Problem happened with multiple resolutions, check the bw file")
  }

  if (length(params_res) == 1){
    params_res = as.integer(params_res)
  } else {
    stop("Problem happened with multiple resolutions (2), check the bw file")
  }

  # now we kill the wrong resolution
  ol = length(gr)
  gr = gr[GenomicRanges::width(gr) == params_res]
  el = length(gr)
  per = scales::percent(el/ol)
  message(glue::glue("After removing not uniform resolution values, {per} regions kept."))


  ### STEP1
  # Should be 1 if it is signal/input, then it is equivalent to log(1) = 0
  # so if the score in in log(fc) this value should be 0.
  params_min_value = min_value

  # these are the bins which don't account for the minimum valkue, thus
  # they are sent to the bin 0! (aka baseline)
  dat_bin0 = gr[gr$score <= params_min_value]

  gr_filt = gr[gr$score > params_min_value]
  # gr_filt$score %>% log %>% hist

  ### STEP2
  # we order the bins according to their score
  gr_filt = gr_filt[order(gr_filt$score)]

  ### STEP3
  # then we bin the genome according to the argument number of bins/
  # because the resolution is shared they will have the same genome size.
  labels = glue::glue("eqFreqBin{seq_len(n_bins)}of{n_bins}")
  bin_values = cut(seq_len(length(gr_filt)),
                   breaks = n_bins,
                   labels = labels)

  gr_filt$bins = bin_values

  ## Now I introduce the info in the 0 bins
  dat_bin0$bins = glue::glue("eqFreqBin0of{n_bins}")

  dat_result = c(dat_bin0,gr_filt)
  dat_result = GenomicRanges::sort(dat_result)
  dat_result
}




binSignal_cumSum <- function(gr,
                             min_value = 1,
                             n_bins = 4) {

  ## see comments above

  params_min_value = min_value

  dat_bin0 = gr[gr$score <= params_min_value]
  gr_filt = gr[gr$score > params_min_value]

  gr_filt = gr_filt[order(gr_filt$score)]
  w_values = GenomicRanges::width(gr_filt)
  gr_filt$cumWidth = cumsum(as.numeric(w_values))

  max_genone = gr_filt[length(gr_filt)]$cumWidth

  # this generates equal sized bins from the minimum value to the maximum
  cut_points = seq(from = params_min_value,
                   to = max_genone,
                   ## I think this + 1 is because these are cutpoints
                   ## so it will generate one bin less.
                   length.out = n_bins + 1)

  labels = glue::glue("eqFreqBin{seq_len(n_bins)}of{n_bins}")
  bin_values = cut(gr_filt$cumWidth,
                   breaks = n_bins,
                   labels = labels)

  gr_filt$bin_values = bin_values

  ## dat_bin0 is a gr object
  if (length(dat_bin0) == 0){
    dat_result = gr_filt
  } else {
    dat_bin0$bin_values = glue::glue("eqFreqBin0of{n_bins}")
    dat_result = c(dat_bin0,gr_filt)
  }

  dat_result = GenomicRanges::sort(dat_result)
  dat_result

}



parallel_intersect <- function(gr_score,gr_mask){

  ## this functions takes a gr object with an score column
  ## (coming from a bw file)
  ## and removes the parts which are not in the
  ## gr_mask by parallel intersection
  ## Keeping original scores to the split ranges.

  if (!any(seqlevels(gr_score) %in% seqlevels(gr_mask))){
    stop("No shared seqlevels")
  }

  mcols(gr_mask) = NULL

  if (!isDisjoint(gr_mask)){
    stop("Mask should not have overlaping segments")
  }

  if (!isDisjoint(gr_score)){
    stop("GR score should not have overlaping segments")
  }

  ## both mask and score gr objects should be reduced here.
  ## Mask we can do but score?

  ovr = findOverlaps(query = gr_score,
                     subject = gr_mask)

  ol = sum(GenomicRanges::width(gr_score))
  # I am pretty sure duplicated overlaps are not a problem, how to test.
  gr_score_int = GenomicRanges::pintersect(x = gr_score[queryHits(ovr)],
                                           y = gr_mask[subjectHits(ovr)])
  el = sum(GenomicRanges::width(gr_score_int))

  if(options()[["verbose"]]){
    message(glue::glue("After the pintersect, {scales::percetn(el/ol)} size remaining"))
  }

  # I thought I had to assign the score of the whole initial bin
  # to the resulting intersected bin but it does so automatically
  # maybe could drop this in function test and move it in the
  # test sections
  score_original = gr_score$score
  resulting_score = score_original[queryHits(ovr)]
  stopifnot(all(resulting_score == gr_score_int$score))

  return(gr_score_int)

}

