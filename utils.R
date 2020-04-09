library(tidyverse)
library(memoise)

set.seed(5)
global_n_grid <- 11
global_n_trials <- 1e3
global_n_patients <- c(12, 24, 48, 96, 192)
global_n_donors <- c(2, 4, 6, 12)
