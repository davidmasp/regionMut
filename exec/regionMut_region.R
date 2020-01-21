#!/usr/bin/env Rscript

# Name: regionMut region
# Author: DMP
# Description: Perform interactions of bed files and count ntps at risk

# params ------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-b", "--regions"),
              action="store",
              type='character',
              help="regions table input (bed files) [default %default]"),
  make_option(c("-k", "--kmer"),
              action="store",
              default = 1,
              type='integer',
              help="kmer of the mutation type [default %default]"),
  make_option(c("-g", "--genome"),
              action="store",
              default = "Hsapiens.UCSC.hg19",
              type='character',
              help="genome alias installed in ur system [default %default]"),
  make_option(c("-r", "--mutRef"),
              action="store",
              default = "C,A",
              type='character',
              help="Comma separated mutation reference set [default %default]"),
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
  make_option(c("-v", "--verbose"),
              action="store_true",
              default=TRUE,
              help="verbosity [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if (interactive()){
  opt$regions = "inst/testdata/channels_bins.tsv"
  opt$regions = "rpTime_rpStrand.tsv"
  opt$folder = "inst/testdata/temp_files"
}

# imports -----------------------------------------------------------------

suppressPackageStartupMessages(library(GenomicRanges))
library(magrittr)
library(regionMut)

# data --------------------------------------------------------------------

col_types = readr::cols(
  master_name = readr::col_character(),
  bins = readr::col_character(),
  file_name = readr::col_character(),
  strand = readr::col_character()
)
input_ref_data = readr::read_tsv(opt$regions,
                           col_types = col_types)


lor_input = read_input_regions(ref_data = input_ref_data)

# script ------------------------------------------------------------------

regions = process_region_interactions2(list_of_regions = lor_input)
genome = helperMut::genome_selector(alias = opt$genome)

## DROP UNUSED REGIONS?
mask = lengths(regions) == 0
warning(glue::glue("{sum(mask)} interactions from the {length(mask)} possible dropped due to lack of interaction"))
regions[mask] = NULL

names(regions) = purrr::map_chr(regions, function(x){
  unique(as.character(x$id))
})

# first we don't simplify the set because we may have strandness
offset_df = get_offset(gr_list = regions,
                       k = opt$kmer,
                       genome = genome,
                       refSet = Biostrings::DNA_BASES)

ref_set = unlist(strsplit(opt$mutRef,split = ","))
stopifnot(length(ref_set) == 2)


## extract metadata from regions


std_tmp = purrr::map_chr(regions,function(x){
  as.character(unique(strand(x)))
  })
strandLess = all(unique(std_tmp) == "*") & (length(unique(std_tmp)) == 1)

## extract and join metadata which is stored in regions
purrr::map_df(regions,function(x){
  dplyr::distinct(as.data.frame(mcols(x)))
}) -> r_ids

r_ids$id  = as.character( r_ids$id )
dplyr::full_join(x = offset_df,
                 y =  r_ids,
                 by = c("combination_id" = "id")) -> offset_df_all


if (strandLess){
  group_vars = unique(input_ref_data$master_name)
  offset_df_all$ctx_simplified = helperMut::simplify_ctx(
                                        as.character(offset_df_all$ctx),
                                        simplify_set = ref_set)
  offset_df_all %>%
    dplyr::group_by_at(c(group_vars,"ctx_simplified")) %>%
    dplyr::summarise(ctx_counts_all = sum(ctx_counts)) -> offset_df_all_simp
} else {
  ## This is only if strand is enforced
  strand_mask = input_ref_data$strand == "*"
  group_vars = unique(input_ref_data[strand_mask,][["master_name"]])
  group_vars_std = unique(input_ref_data[!strand_mask,][["master_name"]])
  message("Strand mode in")
  stopifnot(length(group_vars_std) == 1)

  offset_df_all[[group_vars_std]] %>%
    stringr::str_extract("(?<=_)[:alnum:]+") -> std_group

  offset_df_all[[group_vars_std]] %>%
    stringr::str_extract("[:alnum:]+(?=_)") -> feat_group

  feat_group = unique(feat_group)
  stopifnot(length(feat_group) == 1)

  offset_df_all$ctx_simplified =
    helperMut::simplify_ctx(ctx = as.character(offset_df_all$ctx),
                            simplify_set = ref_set)


  feat_name_true = glue::glue("{group_vars_std}_{feat_group}_ref")
  feat_name_false = glue::glue("{group_vars_std}_anti{feat_group}_ref")
  K = unique(nchar(as.character(offset_df_all$ctx)))
  k = (K-1)/2

  offset_df_all[[group_vars_std]] = ifelse(
    substr(x = offset_df_all$ctx, k + 1, k + 1) %in% ref_set,
    feat_name_true,
    feat_name_false
  )

  offset_df_all %>% dplyr::ungroup() %>%
    dplyr::group_by_at(c(group_vars,group_vars_std,"ctx_simplified")) %>%
    dplyr::summarise(ctx_counts_all = sum(ctx_counts)) -> offset_df_all_simp
}

# output ------------------------------------------------------------------

## make sure folder exists
fs::dir_create(opt$folder)

saveRDS(object = regions,
        file = fs::path(opt$folder,glue::glue("{opt$prefix}_int_regions.rds")))

opath = fs::path(opt$folder,
                 glue::glue("{opt$prefix}_offset.tsv"))
readr::write_tsv(x = offset_df_all_simp,
                 path = opath)
