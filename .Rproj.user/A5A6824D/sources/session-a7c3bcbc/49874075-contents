# Analysis helpers for plot_coding_mixture.R. Base R only; never refits models.
cm_classes <- c("additive", "recessive", "dominant")
cm_short <- c("add", "rec", "dom")
cm_calls <- c(cm_classes, "ambiguous", "undetected")

cm_scenario <- function(counts) {
  labels <- c("Additive", "Recessive", "Dominant")[counts > 0]
  if (!length(labels)) stop("No generating causal effects.")
  if (length(labels) == 1) paste(labels, "only") else paste(labels, collapse = " + ")
}

cm_coding_map <- function(x) {
  p <- length(x$susie_pip)
  m <- length(x$susie_mix_pip)
  map <- x$mix_to_add
  if (!p || !m || length(map) != m || anyNA(map) ||
      any(map < 1 | map > p | map != floor(map))) stop("Invalid mix_to_add map.")
  coding <- x$mix_coding
  origin <- "saved_labels"
  if (is.null(coding) && !is.null(names(x$susie_mix_pip)) &&
      all(grepl("__(additive|recessive|dominant)$", names(x$susie_mix_pip)))) {
    coding <- sub(".*__", "", names(x$susie_mix_pip))
    origin <- "predictor_names"
  }
  if (is.null(coding)) {
    # Compact saves omit labels but preserve A, R, D column order. Within
    # each block, biological-SNP indices are strictly increasing. A reset
    # after the additive block uniquely identifies the R/D boundary.
    if (m < p || !identical(as.integer(map[seq_len(p)]), seq_len(p))) {
      stop("Cannot reconstruct coding: additive block is not 1:p.")
    }
    tail <- map[-seq_len(p)]
    q <- length(tail)
    resets <- which(diff(tail) <= 0)
    if (length(resets) > 1L) stop("Invalid compact coding-block order.")
    lower <- max(0L, q - p)
    upper <- min(p, q)
    if (length(resets)) lower <- upper <- resets[1]
    # Generating predictor indices provide additional boundary constraints.
    for (k in seq_along(x$true_pos_mix)) {
      pos <- x$true_pos_mix[k] - p
      cl <- x$causal_coding[k]
      if (cl == "recessive") lower <- max(lower, pos)
      if (cl == "dominant") upper <- min(upper, pos - 1L)
      if (cl == "additive" && pos > 0) stop("Additive truth lies outside additive block.")
    }
    if (lower > upper) stop("Coding blocks disagree with generating truth.")
    if (lower != upper) {
      stop("Ambiguous compact coding map: save mix_coding in future simulations.")
    }
    coding <- rep(cm_classes, c(p, lower, q - lower))
    origin <- "reconstructed_blocks"
  }
  if (length(coding) != m || anyNA(coding) || any(!coding %in% cm_classes) ||
      anyDuplicated(map + (match(coding, cm_classes) - 1L) * p)) stop("Invalid coding labels or duplicate SNP/coding pairs.")
  if (!identical(as.character(coding[x$true_pos_mix]), as.character(x$causal_coding)) ||
      !identical(as.integer(map[x$true_pos_mix]), as.integer(x$true_pos))) {
    stop("Predictor coding/map disagrees with generating truth.")
  }
  list(coding = coding, origin = origin)
}

cm_assess_replicate <- function(x, snp_cutoff = .9,
                              thresholds = c(0, .1, .5, .8, .9, .95, .99, 1),
                              n_bins = 10L, tie_tolerance = 1e-12) {
  s <- x$settings
  truth <- x$true_pos_mix
  K <- length(x$true_pos)
  p <- length(x$susie_pip)
  pip <- x$susie_mix_pip
  snp_pip <- x$susie_mix_pip_snp
  if (!K || is.null(s) || length(truth) != K || length(x$causal_coding) != K ||
      anyDuplicated(x$true_pos) || anyDuplicated(truth) ||
      anyNA(c(truth, x$true_pos)) || any(truth < 1 | truth > length(pip) | truth != floor(truth)) ||
      any(x$true_pos < 1 | x$true_pos > p | x$true_pos != floor(x$true_pos)) ||
      anyNA(x$causal_coding) || any(!x$causal_coding %in% cm_classes)) stop("Invalid generating truth.")
  tolerance <- 1e-10
  if (!is.numeric(pip) || !length(pip) || any(!is.finite(pip)) || any(pip < -tolerance | pip > 1 + tolerance) ||
      !is.numeric(snp_pip) || length(snp_pip) != p || any(!is.finite(snp_pip)) ||
      any(snp_pip < -tolerance | snp_pip > 1 + tolerance)) stop("Invalid coding or biological-SNP PIPs.")
  # Grouping fitted component probabilities can produce 1 + machine epsilon.
  # Clip only roundoff-sized excursions, and record them in the audit.
  n_clamped <- sum(pip < 0 | pip > 1) + sum(snp_pip < 0 | snp_pip > 1)
  pip <- pmin(1, pmax(0, pip))
  snp_pip <- pmin(1, pmax(0, snp_pip))
  labels <- cm_coding_map(x)
  coding <- labels$coding
  code <- match(coding, cm_classes)
  positive <- seq_along(pip) %in% truth
  at_causal_snp <- x$mix_to_add %in% x$true_pos
  mass <- vapply(cm_classes, function(cl) sum(pip[coding == cl]), numeric(1))
  n_true <- tabulate(match(x$causal_coding, cm_classes), nbins = 3L)
  # Classification is evaluated at the KNOWN causal SNP. It is conditional
  # on the correct location and is not an end-to-end discovery accuracy.
  confusion <- matrix(0L, nrow = 3L, ncol = 5L,
                      dimnames = list(cm_classes, cm_calls))
  for (k in seq_len(K)) {
    ids <- which(x$mix_to_add == x$true_pos[k])
    scores <- setNames(numeric(3), cm_classes)
    scores[coding[ids]] <- pip[ids]
    if (snp_pip[x$true_pos[k]] < snp_cutoff || max(scores) <= 0) {
      call <- "undetected"
    } else {
      winners <- which(abs(scores - max(scores)) <= tie_tolerance)
      call <- if (length(winners) > 1L) "ambiguous" else cm_classes[winners]
    }
    confusion[x$causal_coding[k], call] <- confusion[x$causal_coding[k], call] + 1L
  }
  # Each predictor's target is 1 only for its exact generating SNP AND coding.
  bin <- pmin(n_bins, floor(pip * n_bins) + 1L)
  bin_id <- (code - 1L) * n_bins + bin
  calibration <- matrix(0, nrow = 3L * n_bins, ncol = 4L,
                        dimnames = list(NULL, c("n", "sum_pip", "n_true", "sum_brier")))
  calibration[, "n"] <- tabulate(bin_id, nbins = nrow(calibration))
  binned <- rowsum(cbind(pip, as.numeric(positive), (pip - positive)^2), bin_id, reorder = FALSE)
  calibration[as.integer(rownames(binned)), -1L] <- binned
  nt <- length(thresholds)
  discovery <- matrix(0, nrow = 3L * nt, ncol = 4L,
                      dimnames = list(NULL, c("n_calls", "n_correct", "n_wrong_coding", "n_wrong_snp")))
  for (cl in 1:3) {
    ids <- which(code == cl)
    bins <- findInterval(pip[ids], thresholds)
    selected <- list(rep(TRUE, length(ids)), positive[ids],
                     at_causal_snp[ids] & !positive[ids], !at_causal_snp[ids])
    for (j in seq_along(selected)) {
      discovery[(cl - 1L) * nt + seq_len(nt), j] <-
        rev(cumsum(rev(tabulate(bins[selected[[j]]], nbins = nt))))
    }
  }
  any_false <- vapply(1:3, function(cl) {
    as.integer(any(pip[code == cl & !positive] >= snp_cutoff))
  }, integer(1))
  metrics <- c(setNames(mass, paste0("mass_", cm_short)),
               setNames(n_true, paste0("true_", cm_short)),
               setNames(any_false, paste0("any_false_", cm_short)),
               setNames(as.numeric(t(confusion)),
                        unlist(lapply(cm_short, function(cl) paste0("conf_", cl, "_", cm_calls)))))
  list(metrics = metrics, calibration = calibration, discovery = discovery,
       map_origin = labels$origin, n_clamped = n_clamped)
}

cm_read_simulations <- function(config) {
  files <- list.files(config$chunk_dir, config$file_pattern, full.names = TRUE)
  if (!length(files)) stop("No matching simulation checkpoints in ", config$chunk_dir)
  pattern <- paste0("^.*_add([0-9]+)_rec([0-9]+)_dom([0-9]+)_n([0-9]+)_L([0-9]+)",
                    "_pve([^_]+)_seed([^_]+)_reps([0-9]+)_chunk([0-9]+)\\.RData$")
  parts <- regmatches(basename(files), regexec(pattern, basename(files)))
  if (any(lengths(parts) != 10L)) stop("Unrecognized simulation chunk filename.")
  info <- as.data.frame(do.call(rbind, lapply(parts, function(z) as.numeric(z[-1]))))
  names(info) <- c("L_add", "L_rec", "L_dom", "n", "L", "pve", "seed_base", "requested", "chunk")
  info$file <- files
  info <- info[info$n == config$n & info$L == config$fit_L & info$pve %in% config$pve_values, ]
  info <- info[order(-info$requested, -as.numeric(file.info(info$file)$mtime)), ]
  if (!nrow(info)) stop("No checkpoints match the selected n, fitted L, and PVE.")
  seen <- new.env(hash = TRUE, parent = emptyenv())
  counts <- new.env(hash = TRUE, parent = emptyenv())
  rows <- audits <- issues <- list()
  for (f in seq_len(nrow(info))) {
    z <- info[f, ]
    message("Reading ", f, "/", nrow(info), ": ", basename(z$file))
    saved <- new.env(parent = emptyenv())
    load(z$file, envir = saved)
    if (!exists("results", saved, inherits = FALSE) || !is.list(saved$results)) stop("Missing results list.")
    audit <- data.frame(file = basename(z$file), saved = length(saved$results), examined = 0L,
                        included = 0L, errors = 0L, duplicates = 0L, invalid = 0L,
                        nonconverged = 0L, excluded_nonconverged = 0L,
                        reconstructed_maps = 0L, roundoff_clamped_replicates = 0L)
    metadata <- values <- list()
    for (i in seq_len(min(length(saved$results), config$max_reps_per_file))) {
      x <- saved$results[[i]]
      audit$examined <- audit$examined + 1L
      issue <- NULL
      if (!is.null(x$error)) {
        audit$errors <- audit$errors + 1L
        issue <- paste(x$error, collapse = "; ")
      } else {
        assessed <- tryCatch({
          s <- x$settings
          if (length(x$seed) != 1L || !is.finite(x$seed) || is.null(s) ||
              !isTRUE(all.equal(as.numeric(c(s$L_add, s$L_rec, s$L_dom, s$n, s$L, s$pve)),
                               as.numeric(z[c("L_add", "L_rec", "L_dom", "n", "L", "pve")])))) {
            stop("Settings/seed disagree with checkpoint condition.")
          }
          if (!isTRUE(all.equal(c(s$min_maf, s$hwe_thresh, s$min_n_rec), c(.05, 1e-8, 5)))) {
            stop("Unexpected genotype-QC settings; do not pool with the default simulations.")
          }
          expected <- c(s$L_add, s$L_rec, s$L_dom)
          if (isTRUE(s$all_additive)) expected <- c(sum(expected), 0, 0)
          if (!identical(as.integer(table(factor(x$causal_coding, levels = cm_classes))),
                         as.integer(expected))) stop("Causal coding counts disagree with settings.")
          cm_assess_replicate(x, config$snp_cutoff, config$thresholds, config$n_bins)
        }, error = identity)
        if (inherits(assessed, "error")) {
          audit$invalid <- audit$invalid + 1L
          issue <- conditionMessage(assessed)
        } else {
          design <- paste0("add", z$L_add, "_rec", z$L_rec, "_dom", z$L_dom)
          key <- paste(design, isTRUE(x$settings$all_additive), z$pve, z$n, z$L, x$seed, sep = "|")
          if (exists(key, seen, inherits = FALSE)) {
            audit$duplicates <- audit$duplicates + 1L
            next
          }
          converged <- x$metrics$converged[match("SuSiE-mix", x$metrics$method)]
          if (length(converged) != 1L || is.na(converged)) stop("Missing mixed-fit convergence status.")
          assign(key, TRUE, seen)
          if (!converged) {
            audit$nonconverged <- audit$nonconverged + 1L
            if (config$exclude_nonconverged) {
              audit$excluded_nonconverged <- audit$excluded_nonconverged + 1L
              next
            }
          }
          ntrue <- assessed$metrics[paste0("true_", cm_short)]
          scenario <- cm_scenario(ntrue)
          metadata[[length(metadata) + 1L]] <- data.frame(
            seed = x$seed, scenario = scenario, pve = z$pve, K = sum(ntrue),
            configuration = paste0("add", ntrue[1], "_rec", ntrue[2], "_dom", ntrue[3]),
            design = design, all_additive = isTRUE(x$settings$all_additive),
            converged = converged, row.names = NULL)
          values[[length(values) + 1L]] <- assessed$metrics
          group <- paste(scenario, z$pve, sum(ntrue), sep = "|")
          old <- if (exists(group, counts, inherits = FALSE)) get(group, counts) else NULL
          if (is.null(old)) {
            assign(group, list(calibration = assessed$calibration, discovery = assessed$discovery), counts)
          } else {
            old$calibration <- old$calibration + assessed$calibration
            old$discovery <- old$discovery + assessed$discovery
            assign(group, old, counts)
          }
          audit$included <- audit$included + 1L
          audit$reconstructed_maps <- audit$reconstructed_maps + (assessed$map_origin == "reconstructed_blocks")
          audit$roundoff_clamped_replicates <- audit$roundoff_clamped_replicates + (assessed$n_clamped > 0)
        }
      }
      if (!is.null(issue)) issues[[length(issues) + 1L]] <- data.frame(
        file = basename(z$file), replicate = i, seed = if (length(x$seed) == 1L) x$seed else NA_real_,
        issue = issue)
    }
    if (length(values)) rows[[length(rows) + 1L]] <- cbind(do.call(rbind, metadata), do.call(rbind, values))
    audits[[f]] <- audit
    rm(saved, metadata, values)
  }
  audit <- do.call(rbind, audits)
  issues <- if (length(issues)) do.call(rbind, issues) else
    data.frame(file = character(), replicate = integer(), seed = numeric(), issue = character())
  write.csv(audit, file.path(config$output_dir, "file_audit.csv"), row.names = FALSE)
  write.csv(issues, file.path(config$output_dir, "excluded_replicates.csv"), row.names = FALSE)
  if (!length(rows)) stop("No usable simulations; see excluded_replicates.csv.")
  calibration <- discovery <- list()
  for (key in ls(counts)) {
    fields <- strsplit(key, "|", fixed = TRUE)[[1]]
    meta <- data.frame(scenario = fields[1], pve = as.numeric(fields[2]), K = as.integer(fields[3]))
    value <- get(key, counts)
    calibration[[key]] <- cbind(meta[rep(1L, 3L * config$n_bins), ],
      coding = rep(cm_classes, each = config$n_bins), bin = rep(seq_len(config$n_bins), 3L),
      as.data.frame(value$calibration), row.names = NULL)
    discovery[[key]] <- cbind(meta[rep(1L, 3L * length(config$thresholds)), ],
      coding = rep(cm_classes, each = length(config$thresholds)), threshold = rep(config$thresholds, 3L),
      as.data.frame(value$discovery), row.names = NULL)
  }
  list(replicates = do.call(rbind, rows), calibration = do.call(rbind, calibration),
       discovery = do.call(rbind, discovery), audit = audit, issues = issues)
}

cm_ratio <- function(a, b) ifelse(b > 0, a / b, NA_real_)

cm_mixture_summary <- function(replicates, bootstrap_reps = 500L, seed = 20260911L) {
  keys <- c("scenario", "configuration", "design", "all_additive", "pve", "K")
  groups <- split(seq_len(nrow(replicates)), interaction(replicates[keys], drop = TRUE))
  seeds <- sort(unique(replicates$seed))
  set.seed(seed)
  # Seeds are reused across configurations and PVE values: share the same
  # bootstrap weights instead of treating these rows as independent draws.
  weights <- rmultinom(bootstrap_reps, length(seeds), rep(1, length(seeds)))
  result <- lapply(groups, function(ids) {
    d <- replicates[ids, ]
    mass <- as.matrix(d[paste0("mass_", cm_short)])
    by_seed <- rowsum(mass, d$seed)
    totals <- matrix(0, length(seeds), 3L)
    totals[match(rownames(by_seed), as.character(seeds)), ] <- by_seed
    boot <- crossprod(weights, totals)
    boot <- boot / rowSums(boot)
    estimate <- colSums(mass) / sum(mass)
    truth <- colSums(d[paste0("true_", cm_short)])
    truth <- truth / sum(truth)
    bounds <- vapply(1:3, function(j) {
      valid <- boot[is.finite(boot[, j]), j]
      if (length(unique(d$seed)) < 2L || length(valid) < 2L) return(c(NA_real_, NA_real_))
      unname(quantile(valid, c(.025, .975)))
    }, numeric(2))
    data.frame(d[rep(1L, 3L), keys], coding = cm_classes,
               true_share = unname(truth), estimated_share = unname(estimate),
               bias = unname(estimate - truth), lower = bounds[1, ], upper = bounds[2, ],
               pip_mass = colSums(mass), n_replicates = nrow(d), n_seed_blocks = length(unique(d$seed)),
               zero_pip_replicates = sum(rowSums(mass) == 0), row.names = NULL)
  })
  do.call(rbind, result)
}

cm_confusion_summary <- function(replicates) {
  keys <- c("scenario", "pve", "K")
  vars <- unlist(lapply(cm_short, function(cl) paste0("conf_", cl, "_", cm_calls)))
  pooled <- aggregate(replicates[vars], replicates[keys], sum)
  rows <- lapply(seq_len(nrow(pooled)), function(i) {
    data.frame(pooled[rep(i, 15L), keys], true_coding = rep(cm_classes, each = 5L),
               called_coding = rep(cm_calls, 3L), count = as.numeric(pooled[i, vars]), row.names = NULL)
  })
  do.call(rbind, rows)
}

cm_calibration_rates <- function(table, groups = c("scenario", "pve", "K", "coding", "bin")) {
  out <- aggregate(table[c("n", "sum_pip", "n_true", "sum_brier")], table[groups], sum)
  out$mean_pip <- cm_ratio(out$sum_pip, out$n)
  out$observed_frequency <- cm_ratio(out$n_true, out$n)
  out$brier <- cm_ratio(out$sum_brier, out$n)
  out
}

cm_discovery_rates <- function(table, replicates) {
  out <- table
  groups <- c("scenario", "pve", "K")
  truth <- aggregate(replicates[paste0("true_", cm_short)], replicates[groups], sum)
  key <- function(d) do.call(paste, c(d[groups], sep = "|"))
  row <- match(key(out), key(truth))
  out$n_true <- as.matrix(truth[paste0("true_", cm_short)])[cbind(row, match(out$coding, cm_classes))]
  out$precision <- cm_ratio(out$n_correct, out$n_calls)
  out$recall <- cm_ratio(out$n_correct, out$n_true)
  out$fdp <- cm_ratio(out$n_wrong_coding + out$n_wrong_snp, out$n_calls)
  out
}
