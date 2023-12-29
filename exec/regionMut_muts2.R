#!/usr/bin/env Rscript

# Name: regionMut muts (per folder)
# Author: DMP
# Description: Count mutations in the given region sites


# params ------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-m", "--mutations"),
              action="store",
              type='character',
              help="a folder with vcf, vranges or tsv files [default %default]"),
  make_option(c("-t", "--type"),
              action="store",
              type='character',
              help="the type of files provided [default %default]"),
  make_option(c("-r", "--regions"),
              action="store",
              type='character',
              help="list of regions comming from region subcommand [default %default]"),
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
  make_option(c("-R", "--mutRef"),
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
  make_option(c("-S", "--filterSet"),
              action="store",
              default = NULL,
              type='character',
              help="A mutation set to filter the input table. Use format TCW_T. Multiple sets can be included as colon separated argument TCW_T:NAT_T.[%default]"),
  make_option(c("-v", "--verbose"),
              action="store_true",
              default=FALSE,
              help="verbosity [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

options(verbose = opt$verbose)

if(interactive()){
  opt$regions = "LMR_brain_int_regions.rds"
  opt$mutations = "./data/"
  opt$folder = "inst/testdata/temp_files"
  opt$type = "vcfAll"
}

# imports -----------------------------------------------------------------

suppressPackageStartupMessages(library(VariantAnnotation))
library(helperMut)
library(GenomicRanges)
library(magrittr)
library(regionMut)

# data --------------------------------------------------------------------

## get files:
gtype = switch (opt$type,
  "vcfz" = "*.vcf.gz",
  "vcf" = "*.vcf",
  "vcfAll" = c("*.vcf.gz","*.vcf"),
  "vr" = "*.rds"
)

getFn <- function(x) {
  if (endsWith(x, "vcf.gz")){
    sname_res = stringr::str_extract(x, "[:graph:]+(?=.vcf.gz)")
  } else if (endsWith(x, "vcf")) {
    sname_res = stringr::str_extract(x, "[:graph:]+(?=.vcf)")
  } else if (endsWith(x, "rds")) {
    sname_res = stringr::str_extract(x, "[:graph:]+(?=.rds)")
  } else {
    stop("error!")
  }
  return(sname_res)
}

gtype %>% purrr::map(function(x){
  ifiles = fs::dir_ls(opt$mutations, glob = x)
  ifiles
}) %>% unlist() -> files

names(files) = NULL

files %>% purrr::map(function(fn){
  extension = tools::file_ext(fn)
  sample_name = getFn(fs::path_file(fn))
  if (extension == "rds"){
    dat_vr = readRDS(fn)
  } else {
    dat_vcf = readVcf(fn)
    dat_vr = as(dat_vcf,"VRanges")
  }

  sampleNames(dat_vr) %>% as.character() %>% unique() -> suin
  stopifnot(length(suin) == 1)
  sampleNames(dat_vr) = sample_name

  cnames = colnames(mcols(dat_vr))
  if (any(cnames == "MSks")) {
    message("Using MS provided")
    MS_from_VR = mcols(dat_vr)$MSks
  }

  mcols(dat_vr) = NULL

  if (any(cnames == "MSks")) {
    mcols(dat_vr)$MSks = MS_from_VR
  }

  dat_vr
}) %>% VRangesList() %>% unlist() -> datvrl

names(datvrl) = NULL

unique(as.character(sampleNames(datvrl))) %>% length() -> len_sam
stopifnot(len_sam == length(files))

# rm(dat_vr)
dat_vr = datvrl

sample_lvl = sampleNames(dat_vr) %>% levels()
opt$nSamples = length(files)

## indel filter
nchar_ref = nchar(ref(dat_vr))
nchar_alt = nchar(alt(dat_vr))
indel_mask  = ref(dat_vr) == "-" | nchar_alt > 1 | nchar_ref > 1
ol = length(dat_vr)
dat_vr = dat_vr[!indel_mask]
el= length(dat_vr)
warning(glue::glue("Indel positions removed, {scales::percent(1 - (el/ol))} positions removed" ))

### just to check
# should check for shitty double sampled vcfs
if(any(grepl(pattern = "NORMAL",x = sampleNames(dat_vr)))){
  stop("You have a normal sample in your vcf")
}

N_samples = 1

regions = readRDS(opt$regions)

# script ------------------------------------------------------------------

## filtering chromosome + names
genome = helperMut::genome_selector(alias = opt$genome)

reg_sqls = seqlevelsStyle(regions[[1]])
muts_sqls = seqlevelsStyle(dat_vr)
genome_sqls = seqlevelsStyle(genome)

if (genome_sqls != reg_sqls){
  stop("wrong genome")
} else if (reg_sqls != muts_sqls){
  warning("Switching Seq Level styles for mutations")
  seqlevelsStyle(dat_vr) = reg_sqls
}

# mutation work -----------------------------------------------------------

regions_grl = GRangesList(regions)
sqlevels = intersect(seqlevels(regions_grl),seqlevels(dat_vr))
ol = length(dat_vr)
dat_vr = dat_vr[seqnames(dat_vr) %in% sqlevels]
el = length(dat_vr)

perc = scales::percent(el/ol)
message(glue::glue("{perc} of mutations kept by seqlevels"))

## assigning mutations
ovr = findOverlaps(query = dat_vr,
                   subject = regions_grl)

## subset the mutations inside the regions
ol = length(dat_vr)
dat_vr = dat_vr[queryHits(ovr)] # [^932172]
el = length(dat_vr)
perc = scales::percent(el/ol)
warning(glue::glue("{perc} of mutations found in regions"))

## also drop unused chr
# I think the problem is when the chr are in the levels of regions but not
# in the dat_vr ones. (long shot)
dat_vr = keepSeqlevels(dat_vr,value = sqlevels)

strand_ref_vec = unlist(lapply(regions_grl, function(x){
  unique(as.character(strand(x)))
}))
strand_reg = strand_ref_vec[subjectHits(ovr)]
## dat_vr is already ordered from [^932172]
strand(dat_vr) = strand_reg

## obtain MS
cnames = colnames(mcols(dat_vr))
print(cnames)
if (any(cnames == "MSks")) {
  message("Using MS provided")
  MS = mcols(dat_vr)$MSks
} else {
  MS = helperMut::get_MS_VR(x = dat_vr,
                            k = opt$kmer,
                            genome = genome,
                            keep_strand = TRUE)
}

purrr::map_df(regions,function(x){
  df = dplyr::distinct(as.data.frame(mcols(x)))
  ## there is an issue with names that have special characters such as +
  ## when they are transformed to data.frame, the character goes to . and
  ## then it becomes problematic
  if (any(colnames(df) != colnames(mcols(x)))){
    colnames(df) = colnames(mcols(x))
  }
  df
}) -> rg_id


# I think it was wrong before.
#
# This test should make sense so the id matches the name of the region
# in the Granges List (for the overlap to happen correctly).
stopifnot(all(names(regions_grl) == as.integer(rg_id$id)))

## when there's no mutations in a comb feature then it will not generate
## the feature lvels which is good, however, it will generate NAs
## downstream [^jhsaldgas]
features_df = rg_id[subjectHits(ovr),]
rownames(features_df) = NULL

features_df$id = factor(features_df$id, levels = rg_id$id)

MS = factor(MS,
            levels = helperMut::generate_mut_types(
              k = opt$kmer,
              simplify_set = Biostrings::DNA_BASES))

samples = as.character(sampleNames(dat_vr)) %>%
  factor(levels = sample_lvl)

# I need to do this because the 0 are not couted otherwise
table(MS, sample = samples, id = features_df$id) %>% as.data.frame() -> MS_df

rg_id$id = factor(rg_id$id)
dplyr::full_join(MS_df,rg_id) -> res_df

## there were NAs coming from  [^jhsaldgas], should be solved but check?
stopifnot(!any(is.na(res_df$MS)))

# prepare refset
ref_set = unlist(strsplit(opt$mutRef,split = ","))
stopifnot(length(ref_set) == 2)

## Now grouping by strand
strandLess = all(strand(dat_vr) == "*")

if (strandLess){
  cnames = colnames(res_df)
  # this is the name of the feature
  cnames = cnames[5:ncol(res_df)]
  group_vars = c("sample",cnames)
  res_df$ms_simplified = helperMut::simplify_muts(as.character(res_df$MS),
                                                  simplify_set = ref_set)
  res_df %>% dplyr::group_by_at(c(group_vars,"ms_simplified")) %>%
    dplyr::summarise(ms_counts_all = sum(Freq)) -> res_df_all_simp
} else{
  stop("NOT IMPLEMENTED!")
}

# this is 1!
res_df_all_simp$N_samples = N_samples

# filter by mutation set --------------------------------------------------

if (!is.null(opt$filterSet)){
  set_input = unlist(strsplit(x = opt$filterSet, split = ":"))
  mutSet = helperMut::make_set(set_input)
  res_df_all_simp %<>% dplyr::filter(ms_simplified %in% mutSet)
}

# split table -------------------------------------------------------------

res_df_all_simp %>% split(.$sample) %>%
  purrr::map(function(df){
    sname = unique(df$sample) %>% as.character()
    df %<>% dplyr::ungroup() %>% dplyr::select(-sample)
    fs::dir_create(opt$folder)
    opath = fs::path(opt$folder,
                     glue::glue("{sname}_{opt$prefix}_counts.tsv"))
    readr::write_tsv(x = df,
                     file = opath)
  })

# output ------------------------------------------------------------------



