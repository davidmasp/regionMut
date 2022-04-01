#!/usr/bin/env Rscript

# Name: regionMut_flanks.R
# Author: DMP
# Description: Calculate the flanks from a given bed file

# imports -----------------------------------------------------------------

require(magrittr)
suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-i", "--bed"),
              action="store",
              type='character',
              help="A bed file with genome-wide regions of interest [default %default]"),
  make_option(c("-m", "--mask"),
              action="store",
              type='character',
              help="a mask bed file (inclusion) [default %default]"),
  make_option(c("-g", "--genome"),
              action="store",
              default = "Hsapiens.UCSC.hg19",
              type='character',
              help="genome alias installed in ur system [default %default]"),
  make_option(c("-o", "--overlaps"),
              default = FALSE,
              action="store_true",
              help = "Allow overlaps in the flanked regions, this means that
              the null region will potentially have less coverage than the
              signal"),
  make_option(c("-u", "--onlyUpstream"),
              default = FALSE,
              action="store_true",
              help = "Generate only upstream flank"),
  make_option(c("-d", "--onlyDownstream"),
              default = FALSE,
              action="store_true",
              help = "Generate only downstream flank"),
  make_option(c("-r", "--reduce"),
              action="store_true",
              default = FALSE,
              help = "If input region is not disjoint, reduce it"),
  make_option(c("-s", "--strand"),
              action="store_true",
              default = FALSE,
              help = "Force output strand."),
  make_option(c("-p", "--prefix"),
              action="store",
              default = "output",
              help = "Prefix for regions."),
  make_option(c("-v", "--verbose"),
              default = FALSE,
              action="store_true",
              help = "If input region is not disjoint, reduce it")
)

opt = parse_args(OptionParser(option_list=option_list))

# params ------------------------------------------------------------------

if (interactive()){
  warning("interactive mode!")
  opt$mask = "test_mask.bed"
  opt$bed = "promoters.bed"
  opt$prefix = "outDown"
  opt$reduce = TRUE
  opt$overlaps = TRUE
  opt$onlyUpstream = FALSE
  opt$onlyDownstream = TRUE
}

# data --------------------------------------------------------------------

dat_gr = rtracklayer::import(opt$bed)
mask_gr = rtracklayer::import(opt$mask)

# script ------------------------------------------------------------------


## check if input bed is reduced?
if ( (!GenomicRanges::isDisjoint(dat_gr)) & (!opt$reduce) ){
  stop("Error: Input file is not reduced, use option -r")
} else if ((!GenomicRanges::isDisjoint(dat_gr)) & (opt$reduce) ){
  warning("Reducing overlaps in input file")
  dat_gr = GenomicRanges::reduce(dat_gr)
}

stopifnot(GenomicRanges::isDisjoint(dat_gr))

## remove unmasked items (inclusion)
## currently mask is strandless
stopifnot(GenomicRanges::isDisjoint(mask_gr))
ovr = GenomicRanges::findOverlaps(dat_gr, mask_gr, ignore.strand=TRUE)
ol = length(dat_gr)
dat_gr = dat_gr[S4Vectors::queryHits(ovr)]
el = length(dat_gr)
if (opt$verbose & el != ol) {
  prct = scales::percent(el / ol )
  mssg = glue::glue("Input file filtered by mask: {prct}")
  warning(mssg)
}

width_vec = GenomicRanges::width(dat_gr)
extension_sides = round(width_vec/2)
down_ext = GenomicRanges::resize(dat_gr,
                                 width = width_vec*2,
                                 fix = "start")
up_ext = GenomicRanges::resize(down_ext,
                               width = width_vec*3,
                               fix = "end")

if (!opt$onlyUpstream & !opt$onlyDownstream) {
  ## default option
  flank_extension = extension_sides
  flank_up = GenomicRanges::flank(up_ext, flank_extension,  start = TRUE)
  flank_down = GenomicRanges::flank(up_ext, flank_extension, start = FALSE)
} else if (opt$onlyUpstream & !opt$onlyDownstream) {
  flank_extension = width_vec
  flank_up = GenomicRanges::flank(up_ext, flank_extension,  start = TRUE)
  flank_down = GenomicRanges::GRanges()
} else if (!opt$onlyUpstream & opt$onlyDownstream) {
  flank_extension = width_vec
  flank_up = GenomicRanges::GRanges()
  flank_down = GenomicRanges::flank(up_ext, flank_extension, start = FALSE)
} else {
  stop("it seems that both onlyUpstream and onlyDownstream are set,
       this is not possible.")
}

null_gr = c(flank_up, flank_down)
null_gr = GenomicRanges::sort(null_gr)

if (opt$strand) {
  stop("this is not currently implemented")
} else {
  GenomicRanges::strand(null_gr) = "*"
  GenomicRanges::strand(dat_gr) = "*"
  ## the output is strandless and reduced
  null_gr = GenomicRanges::reduce(null_gr)
  dat_gr = GenomicRanges::reduce(dat_gr)
}

if (opt$overlaps) {
  null_gr2 = GenomicRanges::setdiff(null_gr, dat_gr ,ignore.strand=TRUE)
  null_gr3 = GenomicRanges::intersect(null_gr2, mask_gr,ignore.strand=TRUE)
} else {
  null_gr3 = null_gr
}

stopifnot(GenomicRanges::isDisjoint(c(null_gr3, dat_gr)))

width = GenomicRanges::width
if (opt$verbose){
  pct = scales::percent(sum(width(null_gr3)) / sum(width(dat_gr)))
}

# output ------------------------------------------------------------------

output_fn_signal = glue::glue("{opt$prefix}_signal.bed")
dat_gr %>% rtracklayer::export(output_fn_signal)
output_fn_flanks = glue::glue("{opt$prefix}_flanks.bed")
null_gr3 %>% rtracklayer::export(output_fn_flanks)

