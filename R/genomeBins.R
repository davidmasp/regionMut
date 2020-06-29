# Bining






# ALG)
#
#   0) pre-processing
#   1) Separate the 0 values
#   2) sort
#   3) make a quantiles using indeces
#

# data --------------------------------------------------------------------
library(magrittr)
library(GenomicRanges)
fname = "~/../Desktop/wgEncodeBroadHistoneNhekCtcfStdSig.bigWig"
dat_bw <- rtracklayer::import.bw(fname)


# solve resolution prob ---------------------------------------------------
params_resolutions = table(width(dat_bw))

params_minInstances = 10

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
ol = length(dat_bw)
dat_bw = dat_bw[width(dat_bw) == params_res]
el = length(dat_bw)
per = scales::percent(el/ol)
message(glue::glue("After removing not uniform resolution values, {per} samples kept."))


### STEP1
# this is equivalent to log(1) = 0
# I think this should be like signal / input so 1 is the new 0
params_min_value = 1

dat_bin0 = dat_bw[dat_bw$score <= params_min_value]

dat_bw_filt = dat_bw[dat_bw$score > params_min_value]
dat_bw_filt$score %>% log %>% hist


### STEP2
dat_bw_filt = dat_bw_filt[order(dat_bw_filt$score)]

### STEP3
n_bins = 4
labels = glue::glue("eqFreqBin{seq_len(n_bins)}of{n_bins}")
bin_values = cut(seq_len(length(dat_bw_filt)),breaks = n_bins,labels = labels)
dat_bw_filt$bins = bin_values

dat_bw_filt$width = width(dat_bw_filt)

library(ggplot2)
dat_bw_filt %>% as.data.frame() %>%
  ggplot(aes(x = score,
             fill = bin_values)) +
  geom_histogram(bins = 100) +
  scale_x_log10()

dat_bw_filt_or %>% as.data.frame() %>%
  ggplot(aes(x = score,
             fill = bins)) +
  geom_histogram(bins = 100) +
  scale_x_log10()


# first step here is the intersect

params_min_value = 1
dat_bin0 = dat_bw[dat_bw$score <= params_min_value]
dat_bw_filt = dat_bw[dat_bw$score > params_min_value]

dat_bw_filt = dat_bw_filt[order(dat_bw_filt$score)]
dat_bw_filt$cumWidth = cumsum(width(dat_bw_filt))

max_genone = dat_bw_filt[length(dat_bw_filt)]$cumWidth

n_bins = 4
cut_points = seq(from = params_min_value,
                 to = max_genone,
                 length.out = n_bins + 1)
labels = glue::glue("eqFreqBin{seq_len(n_bins)}of{n_bins}")
bin_values = cut(dat_bw_filt$cumWidth,breaks = n_bins,labels = labels)
table(bin_values)
dat_bw_filt$bin_values = bin_values


dat_bw_filt_or = dat_bw_filt

# From bed to bw
#
# I didn't find a way to do this in bedtools but I may have missed it.
#
#   1) make resolution bins of the genome
#   2) intersect with cds / genes / exons
#   3) report percentatges
#

params_res = 1000
genome = helperMut::genome_selector()
chunks = helperMut::get_region_chunks(gr = genome,wl = params_res)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genes(txdb) -> genes_gr




### make p intersect function
