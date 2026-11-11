# Run these sections separately in RStudio (Ctrl+Enter).
library(readr)
library(ggplot2)

# Paths ----
project <- "/project2/mstephens/wdenault/susie_mix"
root <- file.path(project, "results_slide_prior_em")
iteration <- "iteration_016"
summary_iteration <- "iteration_012"

# Count saved gene results ----
path <- file.path(project, "results_em/iteration_018/results")
lf <- list.files(path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
length(lf)

path <- file.path(root, iteration, "results")
lf <- list.files(path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
length(lf)

# Prior history ----
prior_history <- read_csv(file.path(root, iteration, "prior_history.csv"),
                          show_col_types = FALSE)

# Exclude repeated starting snapshots, but keep the initial online prior.
keep <- is.na(prior_history$stage) | prior_history$stage == 0 |
        prior_history$prior_role == "after_update"
prior_history <- prior_history[keep, ]

tissues <- unique(prior_history$tissue)
weight_columns <- sprintf("w_%02d", 1:17)
slider <- seq(-1, 1, length.out = 17)

par(mfrow = c(5, 6), mar = c(3.5, 4, 2, 1), mgp = c(2.3, 0.65, 0))

for (tissue in tissues) {
  h <- prior_history[prior_history$tissue == tissue, ]
  h <- h[order(h$iteration, h$batch), ]
  w <- as.matrix(h[, weight_columns])
  colours <- hcl.colors(nrow(h), "Dark 3")

  plot(slider, w[1, ], type = "n", ylim = range(w), main = tissue,
       ylab = "Slider prior probability", xlab = "Slider value")
  abline(v = 0, col = "grey80")

  for (k in seq_len(nrow(h))) {
    # Dashed for online updates; solid for full EM and the initial prior.
    online <- h$phase[k] == "minibatch" & h$prior_role[k] == "after_update"
    lines(slider, w[k, ], col = colours[k], lty = if (online) 2 else 1)
  }

  # New histories record whether each online pass contains 9 or 12 batches.
  if (!"batches_in_iteration" %in% names(h)) h$batches_in_iteration <- 9L
  last_batch <- which(h$phase == "minibatch" & h$batch == h$batches_in_iteration &
                      h$prior_role == "after_update")
  for (k in last_batch) {
    lines(slider, w[k, ], col = "black", lwd = 3)
  }
}
par(mfrow = c(1, 1))

# Change between successive priors, separately for each tissue ----
df <- data.frame()

for (tissue in tissues) {
  h <- prior_history[prior_history$tissue == tissue, ]
  h <- h[order(h$iteration, h$batch), ]
  w <- as.matrix(h[, weight_columns])

  for (k in seq_len(nrow(h) - 1L)) {
    df <- rbind(df, data.frame(
      diff = sum(abs(w[k + 1, ] - w[k, ])), tissue = tissue,
      update = k, iteration = h$iteration[k + 1], batch = h$batch[k + 1]
    ))
  }
}

df <- df[order(df$update, match(df$tissue, tissues)), ]
plot(df$diff, ylab = "Sum of absolute changes", xlab = "Comparison index")

ggplot(df, aes(x = update, y = diff, colour = tissue)) +
  geom_line() +
  labs(x = "Prior update", y = "Sum of absolute changes")

# Credible-set delta summaries ----
# FALSE: load your existing summary quickly. TRUE: reread all gene results.
# Rebuilding requires results with alpha_delta still stored in each fit.
rebuild_summary <- FALSE
summary_file <- file.path(root, summary_iteration, "delta_summary.rds")

if (!rebuild_summary) {
  delta_summary <- readRDS(summary_file)
  results <- delta_summary$results
  if (NROW(delta_summary$errors) > 0) {
    print(delta_summary$errors)
    stop("The saved summary contains file errors; inspect them before plotting.")
  }
} else {
  path <- file.path(root, summary_iteration, "results")
  lf <- list.files(path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  if (!length(lf)) stop("No gene results found in ", path)

  # Your original small helper: probabilities conditional on the SNP or CS.
  summarize_delta <- function(weights, grid, prefix) {
    prob <- weights / sum(weights)
    ans <- c(delta_mean = sum(prob * grid),
             p_recessive_like = sum(prob[grid < 0]),
             p_recessive = sum(prob[grid == -1]),
             p_additive = sum(prob[grid == 0]),
             p_dominant_like = sum(prob[grid > 0]))
    setNames(ans, paste0(prefix, names(ans)))
  }

  rows <- list()
  for (k in seq_along(lf)) {
    message(k, "/", length(lf), ": ", basename(lf[k]))
    out <- readRDS(lf[k])
    if (length(out) <= 2L) next  # Keep your original file-level filter.

    for (o in seq_along(out)) {
      if (!is.list(out[[o]])) next
      fit <- out[[o]]$fit_slide_prior
      if (is.null(fit$sets$cs) || !length(fit$sets$cs)) next
      if (is.null(fit$alpha_delta)) {
        stop("alpha_delta is missing in ", basename(lf[k]),
             "; use the saved delta_summary.rds for this iteration.")
      }
      tissue <- if (is.null(names(out))) as.character(o) else names(out)[o]

      for (s in seq_along(fit$sets$cs)) {
        cs <- fit$sets$cs[[s]]
        ell <- fit$sets$cs_index[s]
        alpha_cs <- fit$alpha[ell, cs]
        lead <- cs[which.max(alpha_cs)]

        cs_weights <- apply(fit$alpha_delta[ell, cs, , drop = FALSE], 3, sum)
        cs_summary <- summarize_delta(cs_weights, fit$delta_grid, "cs_")
        lead_summary <- summarize_delta(fit$alpha_delta[ell, lead, ],
                                        fit$delta_grid, "lead_")

        rows[[length(rows) + 1L]] <- data.frame(
          file = basename(lf[k]), tissue = tissue,
          cs = names(fit$sets$cs)[s], component = ell, n_snps = length(cs),
          cs_coverage = sum(alpha_cs), lead_index = lead,
          lead_alpha = unname(fit$alpha[ell, lead]),
          lead_delta_forced = unname(fit$delta_forced[lead]),
          lead_n0 = unname(fit$genotype_counts[lead, "0"]),
          lead_n1 = unname(fit$genotype_counts[lead, "1"]),
          lead_n2 = unname(fit$genotype_counts[lead, "2"]),
          min_obs = fit$min_obs,
          cs_forced_mass = sum(alpha_cs * fit$delta_forced[cs]) / sum(alpha_cs),
          as.list(cs_summary), as.list(lead_summary), row.names = NULL
        )
      }
    }
  }
  results <- if (length(rows)) do.call(rbind, rows) else data.frame()
  # No files are overwritten. Save a new summary yourself if wanted.
}

# Inspect and filter the summary ----
if (!nrow(results)) stop("The summary contains no credible sets.")
head(results)
table(results$lead_delta_forced)
mean(results$lead_delta_forced)

# Logical indexing also works when no leads are forced; -which(...) does not.
unforced <- results[!results$lead_delta_forced, ]
if (!nrow(unforced)) stop("There are no unforced lead SNPs to plot.")

hist(unforced$lead_delta_mean, nclass = 200,
     main = "Lead SNPs without forced additive coding",
     xlab = "Posterior mean delta")

# Tissue-specific densities; red = all tissues pooled ----
df_lead <- data.frame(x = unforced$lead_delta_mean, tissue = unforced$tissue)
df_cs <- data.frame(x = unforced$cs_delta_mean, tissue = unforced$tissue)

P1 <- ggplot(df_lead, aes(x = x)) +
  geom_density() +
  geom_density(data = data.frame(x = df_lead$x), colour = "red") +
  facet_wrap(~ tissue) + geom_vline(xintercept = 0) +
  labs(title = "Lead CS SNP delta", x = "Posterior mean delta")

P2 <- ggplot(df_cs, aes(x = x)) +
  geom_density() +
  geom_density(data = data.frame(x = df_cs$x), colour = "red") +
  facet_wrap(~ tissue) + geom_vline(xintercept = 0) +
  labs(title = "CS delta", x = "Posterior mean delta")

P1
P2
gridExtra::grid.arrange(P1, P2, ncol = 2)

# CS probabilities, keeping CSs whose lead is unforced ----
par(mfrow = c(3, 1))
hist(unforced$cs_p_additive, nclass = 100,
     main = "CS: additive", xlab = "Posterior probability of delta = 0")
hist(unforced$cs_p_recessive_like, nclass = 100,
     main = "CS: recessive-like", xlab = "Posterior probability of delta < 0")
hist(unforced$cs_p_dominant_like, nclass = 100,
     main = "CS: dominant-like", xlab = "Posterior probability of delta > 0")
par(mfrow = c(1, 1))

# Lead probabilities ----
hist(unforced$lead_p_additive, nclass = 100, main = "Lead: P(delta = 0)")
# The saved summary has dominant_like (delta > 0), not dominant (delta = 1).
hist(unforced$lead_p_dominant_like, nclass = 100, main = "Lead: P(delta > 0)")
hist(unforced$lead_p_recessive, nclass = 100, main = "Lead: P(delta = -1)")
hist(unforced$lead_p_recessive_like, nclass = 100, main = "Lead: P(delta < 0)")

quantile(unforced$lead_delta_mean)
quantile(results$lead_delta_mean[
  !(results$lead_delta_mean == 0 & results$component == 1)
])

# Same coding as your original regression: 0 for component 1, 1 for others.
unforced$is_later_component <- as.integer(unforced$component != 1)
summary(lm(lead_delta_mean ~ 1, data = unforced))
summary(lm(lead_delta_mean ~ is_later_component, data = unforced))

hist(unforced$cs_delta_mean, nclass = 100,
     main = "Posterior mean delta per CS", xlab = "Posterior mean delta")
delta_pos <- results$cs_delta_mean
delta_lead <- results$lead_delta_mean

# Excluding the same sex-specific tissues as in your original script ----
sex_specific <- c("Vagina", "Uterus", "Breast", "Ovary", "Prostate")
other_tissues <- unforced[!unforced$tissue %in% sex_specific, ]

hist(other_tissues$lead_delta_mean, nclass = 200,
     main = "Unforced lead SNPs, excluding sex-specific tissues",
     xlab = "Posterior mean delta")

q0 <- quantile(unforced$lead_delta_mean, probs = seq(0, 1, length.out = 100))
q1 <- quantile(other_tissues$lead_delta_mean, probs = seq(0, 1, length.out = 100))
plot(q0, q1, main = "QQ plot: posterior mean delta",
     xlab = "Including sex-specific tissues", ylab = "Excluding sex-specific tissues")
abline(a = 0, b = 1)
