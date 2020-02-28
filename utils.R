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
