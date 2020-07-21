#!/usr/bin/env Rscript

# Name: regionMut bin
# Author: DMP
# Description: Transform a bigWig signal in equal-sized bins

# params ------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-s", "--bigWig"),
              action="store",
              type='character',
              help="A bigWig file with genome-wide scores [default %default]"),
  make_option(c("-m", "--mask"),
              action="store",
              type='character',
              help="a mask bed file (inclusion) [default %default]"),
  make_option(c("-g", "--genome"),
              action="store",
              default = "Hsapiens.UCSC.hg19",
              type='character',
              help="genome alias installed in ur system [default %default]"),
  make_option(c("-n", "--numberBins"),
              action="store",
              default = 4,
              type='integer',
              help="Number of bins to output [default %default]"),
  make_option(c("-l", "--minValue"),
              action="store",
              default = 1,
              type="double",
              help="Minimum value to consider in the signal distribution, values smaller to this parameter will be sent to bin 0 [default %default]."),
  make_option(c("-p", "--prefix"),
              action="store",
              default = "output",
              type='character',
              help="Output prefix [default %default]"),
  make_option(c("-f", "--folder"),
              action="store",
              default = ".",
              type='character',
              help="Output folder [default %default]"),
  make_option(c("-P", "--plot"),
              action="store_true",
              default=FALSE,
              help="plot final bins [default %default]"),
  make_option(c("-v", "--verbose"),
              action="store_true",
              default=TRUE,
              help="verbosity [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if (interactive()){
  opt$bigWig = "inst/testdata/wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig"
  opt$mask = "inst/testdata/genes.bed"
}

# imports -----------------------------------------------------------------

suppressPackageStartupMessages(library(GenomicRanges))
library(magrittr)
library(regionMut)

# data --------------------------------------------------------------------
genome = helperMut::genome_selector(alias = opt$genome)

dat_bw = rtracklayer::import.bw(opt$bigWig)
seqlevelsStyle(dat_bw) = seqlevelsStyle(genome)
seqlevels(dat_bw) = seqlevels(genome)
seqlengths(dat_bw) = seqlengths(genome)

if (!is.null(opt$mask)){
  warning("using mask option")
  dat_mask = rtracklayer::import.bed(opt$mask)
  seqlevelsStyle(dat_mask) = seqlevelsStyle(genome)
  seqlevels(dat_mask) = seqlevels(genome)
  seqlengths(dat_mask) = seqlengths(genome)
}


# script ------------------------------------------------------------------

if (!is.null(opt$mask)){
  dat_bw_new = parallel_intersect(dat_bw,dat_mask)
} else {
  dat_bw_new = dat_bw
}


binSignal_cumSum(gr = dat_bw_new,
                 min_value = opt$minValue,
                 n_bins = opt$numberBins) -> result_bins


result_bins_df = result_bins
result_bins_df$width = width(result_bins_df)
as.data.frame(result_bins_df) -> result_bins_df

# output ------------------------------------------------------------------

## make sure folder exists
fs::dir_create(opt$folder)

# score plot
if (opt$plot){
  require(ggplot2)
  ggplot(data = result_bins_df,
         aes(x = score,
             fill = bin_values)) +
    geom_histogram() +
    theme_classic() -> p1

  result_bins_df %>% dplyr::group_by(bin_values) %>%
    dplyr::summarise(total_width = sum(width)) %>%
    ggplot(aes(x = bin_values,
               y = total_width,
               fill = bin_values)) +
    geom_col() +
    theme_classic() -> p2

  result_bins_df %>%
    ggplot(aes(x = score,
               y = cumWidth,
               color = bin_values)) +
    geom_line(size = 3) +
    labs(y = "Cummulative genome-width (bp)",
         x = "score") +
    theme_classic() -> p3

  require(patchwork)

  fp = p1 + (p2 / p3) + plot_layout(guides = "collect")

  oname = glue::glue("{opt$prefix}_binsPlot.pdf")
  fs::dir_create(opt$folder)
  ggsave(filename = oname,path = opt$folder,plot = fp)

}

# save bins ---------------------------------------------------------------

binTypes = result_bins$bin_values
results_split = split(x = result_bins, f = binTypes)

lapply(results_split, reduce) -> results_split_reduced

lapply(results_split_reduced, function(x){
  binName = as.character(unique(x$bin_values))
  oname = glue::glue("{opt$prefix}_{binName}_lowThr{opt$minValue}.bed.gz")
  opath = fs::path(opt$folder,oname)
  mcols(x) = NULL
  x = sort(sortSeqlevels(x))
  rtracklayer::export.bed(object = x,
                          con = opath)
})


