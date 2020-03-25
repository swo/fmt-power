#!/usr/bin/env Rscript

library(tools)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
url <- args[1]
expected_md5 <- args[2]
dest <- args[3]

download.file(url, dest)

observed_md5 <- md5sum(dest)
if (observed_md5 != expected_md5) {
  stop(str_glue("Observed MD5 {observed_md5} != expected MD5 {expected_md5} for file {dest}"))
}
