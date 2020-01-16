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
              help="A mutation set to filter the input table. Use format TCW_T. Multiple sets can be included as comma separated argument TCW_T,NAT_T.[%default]"),
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
  opt$offset = "output_offset.tsv"
  opt$counts = "output_counts.tsv"
  opt$formula = "inst/testdata/test.yml"
  opt$filterSet = "TCW_T,NAT_T"
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

dat = dplyr::inner_join(counts,
                        offset,
                        by = c(group_vars,"ctx" = "ctx_simplified"))

dat$ln_at_risk = log(dat$N_samples * dat$ctx_counts_all)


# filter by mutation set --------------------------------------------------

if (!is.null(opt$filterSet)){
  set_input = unlist(strsplit(x = opt$filterSet, split = ","))
  mutSet = helperMut::make_set(set_input)
  dat %<>% dplyr::filter(dat$ms_simplified %in% mutSet)
}

# formula info --------------------------------------------------------------

opt$formula

form_yml = yaml::read_yaml(opt$formula)

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
  ci_vars = form_yml$ci$variables
  ci_alpha = form_yml$ci$alpha
} else {
  ci_method = "profile"
  ci_vars = group_vars
  ci_alpha = 0.05
}

glm_nb_wrapper(data = dat,
               formula = formula_str,
               ci_method = ci_method,
               ci_parameters = ci_vars,
               alpha = ci_alpha) -> test_coef

# output ------------------------------------------------------------------

fs::dir_create(opt$folder)
opath = fs::path(opt$folder,
                 glue::glue("{opt$prefix}_coef.tsv"))
readr::write_tsv(x = test_coef,path = opath)

