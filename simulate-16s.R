#!/usr/bin/env Rscript

source("utils.R")
library(vegan)

# Load and process ----------------------------------------------------

load_study_data1 <- function(study_name, metadata_clean_f) {
  otu_fn <- str_glue("data/{study_name}_results/{study_name}.otu_table.100.denovo")
  metadata_fn <- str_glue("data/{study_name}_results/{study_name}.metadata.txt")

  metadata <- read_tsv(metadata_fn) %>%
    metadata_clean_f() %>%
    select(sample_id, case_control) %>%
    arrange(sample_id)

  otu_matrix <- read_tsv(otu_fn)

  # Look for common IDs
  sample_ids <- intersect(metadata$sample_id, names(otu_matrix))

  # Throw out IDs in metadata but not in OTU table
  metadata <- metadata %>% filter(sample_id %in% sample_ids)

  otu_matrix <- otu_matrix %>%
    # Keep only columns with sample IDs in the relevant metadata
    select_at(sample_ids) %>%
    as.matrix() %>%
    # Keep only rows (OTUs) with nonzero counts across these samples
    { .[rowSums(.) > 0, ] }

  dissim <- vegdist(t(otu_matrix))

  list(dissim = dissim, metadata = metadata)
}

load_study_data <- memoise(load_study_data1, cache = cache_filesystem("cache/study_data/"))

# Simulate ------------------------------------------------------------

subset_dissim <- function(dissim, idx) {
  dissim %>%
    as.matrix() %>%
    `[`(idx, idx) %>%
    as.dist()
}

simulate_f <- function(dissim, case_idx, control_idx, n, p1, phi = 0.5) {
  # draw donor qualities
  n_good <- rbinom(1, n, phi)
  n_bad <- n - n_good

  # subset dissim matrix
  good_idx <- sample(control_idx, n_good, replace = TRUE)
  bad_idx <- sample(case_idx, n_bad, replace = TRUE)
  Y <- subset_dissim(dissim, c(good_idx, bad_idx))

  # simulate patient outcomes
  donor_p <- c(rep(p1, n_good), rep(1 - p1, n_bad))
  X <- rbinom(n, 1, donor_p)

  # if outcomes are all 0 or all 1, pval = 1.0
  if (length(unique(X)) == 1) return(1.0)

  # ignore complete enumeration warnings
  adonis(Y ~ X)$aov.tab$`Pr(>F)`[1]
}

results_f <- function(study_name, metadata_clean_f) {
  study_data <- load_study_data(study_name, metadata_clean_f)
  metadata <- study_data$metadata
  dissim <- study_data$dissim

  stopifnot(all(nrow(metadata) == dim(dissim)))

  # Which indices in the original data are controls/cases?
  case_idx <- which(metadata$case_control == "case")
  control_idx <- which(metadata$case_control == "control")

  # multi draw
  simulate_trials <- function(n_patients, effect_size) {
    map_dbl(1:global_n_trials, ~ simulate_f(dissim, case_idx, control_idx, n_patients, effect_size))
  }

  f <- memoise(simulate_trials, cache = cache_filesystem(str_glue("cache/{study_name}")))

  crossing(
    n_patients = global_n_patients,
    effect_size = seq(0.5, 1, length.out = global_n_grid)
  ) %>%
    mutate(
      simulation = study_name,
      n_donors = n_patients,
      p_values = pmap(list(n_patients, effect_size), simulate_trials),
      x = map_dbl(p_values, ~ sum(. <= 0.05)),
      n = map_dbl(p_values, length),
      estimate = x / n
    ) %>%
  select(simulation, n_donors, n_patients, effect_size, x, n, estimate)
}

# Data cleaning functions ---------------------------------------------
cdi_clean <- function(df) {
  df %>%
    filter(DiseaseState %in% c("H", "nonCDI")) %>%
    mutate(case_control = recode(DiseaseState, "H" = "control", "nonCDI" = "case"))
}

ob_clean <- function(df) {
  df %>%
    select(-sample_id) %>%
    rename(sample_id = X1) %>%
    mutate(case_control = recode(DiseaseState, "H" = "control", "OB" = "case"))
}

# Run simulations and save results ------------------------------------
cdi_results <- results_f("cdi_schubert", cdi_clean)
ob_results <- results_f("ob_goodrich", ob_clean)

results <- bind_rows(cdi_results, ob_results)
write_tsv(results, "results/16s.tsv")
