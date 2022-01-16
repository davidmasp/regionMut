#!/usr/bin/env Rscript

# Name: regionMut muts
# Author: DMP
# Description: Count mutations in the given region sites


# params ------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-m", "--mutations"),
              action="store",
              type='character',
              help="a vcf, vranges or tsv file [default %default]"),
  make_option(c("-r", "--regions"),
              action="store",
              type='character',
              help="list of regions comming from region subcommand [default %default]"),
  make_option(c("-N", "--nSamples"),
              action="store",
              default = NULL,
              type='integer',
              help="Number of samples, if null, # of samp in the vcf [default %default]"),
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
  opt$regions = "inst/testdata/temp_files/output_int_regions.rds"
  opt$mutations = "inst/testdata/test_vcf.vcf"
  opt$folder = "inst/testdata/temp_files"
}

# imports -----------------------------------------------------------------

suppressPackageStartupMessages(library(VariantAnnotation))
library(helperMut)
library(GenomicRanges)
library(magrittr)
library(regionMut)

# data --------------------------------------------------------------------

## determine format from the extension

extension = tools::file_ext(opt$mutations)
if (extension == "rds"){
  dat_vr = readRDS(opt$mutations)
} else {
  dat_vcf = readVcf(opt$mutations)
  dat_vr = as(dat_vcf,"VRanges")
}

mcols(dat_vr) = NULL

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


if (length(unique(sampleNames(dat_vr))) > 1){
  if (!is.null(opt$nSamples)){
    warning("Input is not unisample, forcing because -N is set")
    sampleNames(dat_vr) = as.character("UNISAMPLE")
  } else {
    stop("Unisample vcf are required, remove sample info and add -N arg")
  }
}

if (is.null(opt$nSamples)){
  N_samples = 1
} else {
  N_samples = opt$nSamples
}

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
scales::percent(el/ol)

## assigning mutations
ovr = findOverlaps(query = dat_vr,subject = regions_grl)

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
MS = helperMut::get_MS_VR(x = dat_vr,
                     k = opt$kmer,
                     genome = genome,
                     keep_strand = TRUE)

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

features_df$id = factor(features_df$id,levels = rg_id$id)

MS = factor(MS,
            levels = helperMut::generate_mut_types(
              k = opt$kmer,
              simplify_set = Biostrings::DNA_BASES))

# I need to do this because the 0 are not couted otherwise
table(MS,id = features_df$id) %>% as.data.frame() -> MS_df


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
  cnames = cnames[4:ncol(res_df)]
  group_vars = cnames
  res_df$ms_simplified = helperMut::simplify_muts(as.character(res_df$MS),
                                              simplify_set = ref_set)
  res_df %>% dplyr::group_by_at(c(group_vars,"ms_simplified")) %>%
    dplyr::summarise(ms_counts_all = sum(Freq)) -> res_df_all_simp
} else{
  # with strand info
  # find which columns have the strand info
  cnames = colnames(res_df)
  cnames = cnames[4:ncol(res_df)]

  ## strand features need to carry minus
  strand_mask = res_df[,cnames] %>% purrr::map(grepl,
                                               pattern = "minus" ) %>%
    purrr::map_lgl(any)
  group_vars = cnames[!strand_mask]
  group_vars_std = cnames[strand_mask]

  stopifnot(length(group_vars_std) == 1)

  res_df[[group_vars_std]] %>%
    stringr::str_extract("(?<=_)[:alnum:]+") -> std_group

  res_df[[group_vars_std]] %>%
    stringr::str_extract("[:alnum:]+(?=_)") -> feat_group

  feat_group = unique(feat_group)
  stopifnot(length(feat_group) == 1)

  res_df$ms_simplified =
    helperMut::simplify_muts(muts = as.character(res_df$MS),
                             simplify_set = ref_set)

  feat_name_true = glue::glue("{group_vars_std}_{feat_group}_ref")
  feat_name_false = glue::glue("{group_vars_std}_anti{feat_group}_ref")
  k = opt$kmer

  res_df[[group_vars_std]] = ifelse(
    substr(x = res_df$MS, k + 1, k + 1) %in% ref_set,
    feat_name_true,
    feat_name_false
  )

  res_df %>% dplyr::ungroup() %>%
    dplyr::group_by_at(c(group_vars,group_vars_std,"ms_simplified")) %>%
    dplyr::summarise(ms_counts_all = sum(Freq)) -> res_df_all_simp
}

res_df_all_simp$N_samples = N_samples


# filter by mutation set --------------------------------------------------

if (!is.null(opt$filterSet)){
  set_input = unlist(strsplit(x = opt$filterSet, split = ":"))
  mutSet = helperMut::make_set(set_input)
  res_df_all_simp %<>% dplyr::filter(ms_simplified %in% mutSet)
}

# output ------------------------------------------------------------------

fs::dir_create(opt$folder)
opath = fs::path(opt$folder,
                 glue::glue("{opt$prefix}_counts.tsv"))
readr::write_tsv(x = res_df_all_simp,
                 path = opath)

