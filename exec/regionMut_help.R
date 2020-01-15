#!/usr/bin/env Rscript

# Name: regionmut help
# Author: DMP
# Description: List modes of regionmut

# imports -----------------------------------------------------------------

library(magrittr)

# script ------------------------------------------------------------------

path = system.file("exec",package = "regionMut")
available_modes = fs::dir_ls(path) %>%
  fs::path_file() %>%
  stringr::str_extract("(?<=regionMut_)[:alnum:]+(?=.R)")

pv = utils::packageVersion("regionMut")
message(strwrap(glue::glue("regionMut ({pv})"), prefix = " ", initial = ""))

message("Available Modes:")
available_modes %>% purrr::walk(message)


message("Run regionmut MODE -h to obtain help for each mode")
