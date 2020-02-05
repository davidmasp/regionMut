#!/usr/bin/env Rscript

# Name: regionMut regression
# Author: DMP
# Description: Merge counts and contexts, then perform a regression

# params ------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-c", "--counts"),
              action="store",
              type='character',
              help="counts input "),
  make_option(c("-o", "--offset"),
              action="store",
              type='character',
              help="offset input "),
  make_option(c("-F", "--formula"),
              action="store",
              type='character',
              help="Formula specifications in yaml format"),
  make_option(c("-S", "--filterSet"),
              action="store",
              default = NULL,
              type='character',
              help="A mutation set to filter the input table. Use format TCW_T. Multiple sets can be included as colon separated argument TCW_T:NAT_T.[%default]"),
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
              default=FALSE,
              help="verbosity [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

options(verbose = opt$verbose)

if(interactive()){
  opt$offset = "inst/testdata/temp_files/output_offset.tsv"
  opt$counts = "inst/testdata/temp_files/output_counts.tsv"
  opt$formula = "inst/testdata/david_strand_formula.yml"
  opt$filterSet = NULL
}

# imports -----------------------------------------------------------------

library(regionMut)
library(magrittr)

# data --------------------------------------------------------------------

counts = readr::read_tsv(file = opt$counts)
offset = readr::read_tsv(file = opt$offset)

# script ------------------------------------------------------------------

K = unique(nchar(offset$ctx_simplified))
counts$ctx = substr(counts$ms_simplified,1,K)

cnames = colnames(counts)
group_vars = cnames[!cnames %in% c("ms_simplified",
                                   "ms_counts_all",
                                   "N_samples",
                                   "ctx")]
counts %>% dplyr::distinct_at(c(group_vars, "ms_simplified")) %>% nrow() -> urow

if(nrow(counts)> urow){
  ### NOTE:
  ### This is for the case that multiple unisample counts have been concatenated
  ### together. If this is the case, then it summarises in a single table.
  counts %>% dplyr::group_by_at(c(group_vars, "ms_simplified","ctx")) %>%
    dplyr::summarise(
      ms_counts_all = sum(ms_counts_all),
      N_samples = sum(N_samples)
    )
}

dat = dplyr::left_join(counts,
                        offset,
                        by = c(group_vars,"ctx" = "ctx_simplified"))

dat$ln_at_risk = log(dat$N_samples * dat$ctx_counts_all)

# remove impossble regions
dat %<>% dplyr::filter(ctx_counts_all != 0)

# filter by mutation set --------------------------------------------------

# I think this won't be necessary now...
if (!is.null(opt$filterSet)){
  set_input = unlist(strsplit(x = opt$filterSet, split = ":"))
  mutSet = helperMut::make_set(set_input)
  dat %<>% dplyr::filter(dat$ms_simplified %in% mutSet)
}

# formula info --------------------------------------------------------------

form_yml = yaml::read_yaml(opt$formula)


## fix plus signs in the levels
## I think this should be enforced outside the package because it can
## bring problems I think at some point
subs_plus_sign <- function(text) {
  stringr::str_replace_all(
    pattern = "(?<![:blank:])[+](?![:blank:])",
    replacement = ".",
    string = text)
}
if (any(grepl(names(form_yml$levels),pattern = "[+]"))){
  form_yml$formula$variables = subs_plus_sign(form_yml$formula$variables )
  names(form_yml$levels) = subs_plus_sign(names(form_yml$levels))
  colnames(dat) = subs_plus_sign(colnames(dat))
}


value = form_yml$formula$value
variables = form_yml$formula$variables
use_offset = form_yml$formula$offset
if (!use_offset){
  formula_str = glue::glue("{value} ~ {variables}")
} else {
  formula_str = glue::glue("{value} ~ {variables} + offset(ln_at_risk)")
}

for (i in names(form_yml$levels)){
  dat[[i]] = forcats::fct_relevel(dat[[i]],form_yml$levels[[i]])
}

# regression --------------------------------------------------------------

if ("ci" %in% names(form_yml)){
  ci_method = form_yml$ci$method
  ci_alpha = form_yml$ci$alpha
} else {
  ci_method = "profile"
  ci_alpha = 0.05
}

glm_nb_wrapper(data = dat,
               formula = formula_str,
               ci_method = ci_method,
               alpha = ci_alpha) -> test_coef

# output ------------------------------------------------------------------

fs::dir_create(opt$folder)
opath = fs::path(opt$folder,
                 glue::glue("{opt$prefix}_coef.tsv"))
readr::write_tsv(x = test_coef,path = opath)

