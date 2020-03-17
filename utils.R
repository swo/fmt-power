this_file <- function() {
  grep("--file", commandArgs(), value = TRUE) %>%
    str_match("--file=./(.+)") %>%
    `[`(, 2)
}

safe_fisher_test <- function(x) {
  try(return(fisher.test(x)), silent = TRUE)
  chisq.test(x)
}

simulate <- function(f, n_iter, ...) {
  args <- list(...)
  control <- c(args, script = this_file(), n_iter = n_iter)

  hash <- digest::digest(control)
  cache_fn <- str_glue("cache/{hash}.rds")

  if (file.exists(cache_fn)) {
    readRDS(cache_fn)
  } else {
    results <- map(1:n_iter, ~ do.call(f, args))
    saveRDS(results, cache_fn)
    results
  }
}

results_f <- function(effect_size_f) {
  crossing(
    n_donors = c(2, 4, 6, 12, 24),
    n_patients = c(12, 24, 48, 96, 192, 384),
  ) %>%
    filter(n_patients > n_donors) %>%
    mutate(
      patients_per_donor = n_patients / n_donors,
      effect_size = map2_dbl(n_donors, patients_per_donor, effect_size_f)
    )
}

plot_f <- function(results) {
  results %>%
    mutate_at(c("n_donors", "n_patients"), factor) %>%
    ggplot(aes(n_donors, n_patients)) +
    geom_tile(aes(fill = effect_size)) +
    geom_text(aes(label = round(effect_size, 2))) +
    scale_fill_gradient(
      name = "Effect size",
      labels = function(x) round(x, 2),
      low = "white", high = "red"
    ) +
    cowplot::theme_half_open()
}
