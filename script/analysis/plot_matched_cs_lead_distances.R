# Compare ALL gene-tissue regions with the same positive number of credible
# sets in SuSiE and the unweighted SuSiE-mix. Requires only base R.
#
# Run on RCC after generate_summary_results.R:
#   Rscript script/analysis/plot_matched_cs_lead_distances.R
# Other locations / CSV inputs:
#   Rscript script/analysis/plot_matched_cs_lead_distances.R --project-dir=/path/to/project
#   Rscript script/analysis/plot_matched_cs_lead_distances.R --summary-file=res_summary.csv --cs-summary-file=res_cs_summary.csv
# Replot a tissue WITHOUT loading fits or repeating the matching:
#   Rscript script/analysis/plot_matched_cs_lead_distances.R --pairs-file=plot/matched_cs_lead_distances/matched_cs_pairs.csv --tissue=Blood
# Source this file to call match_equal_cs_regions() or plot_matched_cs_distances().
#
# Estimand: one observation per paired CS, so a region with k CSs contributes k
# observations. No association, read-count, convergence, or coding filter is
# applied. Those available region metadata are retained for later filtering.
#
# Matching: sort both models' leads by physical GRCh38 position and pair ranks.
# For absolute distance on one chromosome this globally minimizes the SUM of
# distances over all one-to-one assignments. Independent nearest neighbours
# are inappropriate because they can reuse a CS. Among minimum-total solutions,
# sorted pairing also minimizes the sum of squared distances. This tie rule can
# prefer (100,100) bp over (0,200) bp; it does not maximize unchanged SNPs.
# Remaining ties use SNP identity and CS number for reproducibility. Matching
# is purely positional: it does not establish that CSs are the same signal.
#
# Leads/PIPs are taken from the CS-centric summary: its component-alpha lead
# within each CS (PIP fallback when the component is unavailable). SuSiE-mix PIP
# is predictor/coding-specific, NOT a sum across codings and NOT CS coverage.
# Parse positions from lead_snp itself, never from a TSS distance or another
# saved workhorse SNP. Invalid/incomplete regions are audited, never partially
# matched. Zero-distance pairs are retained and plotted in a separate strip.

matched_cs_require <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing: ", paste(missing, collapse = ", "),
         ". Regenerate the summaries with generate_summary_results.R.")
  }
}

matched_cs_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

matched_cs_key <- function(gene, tissue) {
  # Length prefixes avoid collisions if identifiers themselves contain '|'.
  if (!length(gene)) return(character())
  paste0(nchar(gene), ":", gene, "|", nchar(tissue), ":", tissue)
}

matched_cs_parse_snps <- function(ids) {
  ids <- as.character(ids)
  pattern <- "^(chr(?:[0-9]+|X|Y|M|MT))_([0-9]+)_[^_]+_[^_]+_b38(?:_[^_]+)?$"
  valid <- !is.na(ids) & grepl(pattern, ids, perl = TRUE, ignore.case = TRUE)
  chromosome <- rep(NA_character_, length(ids))
  position <- rep(NA_real_, length(ids))
  chromosome[valid] <- toupper(sub("^chr([^_]+)_.*$", "\\1", ids[valid], ignore.case = TRUE))
  chromosome[chromosome == "MT" & !is.na(chromosome)] <- "M"
  position[valid] <- matched_cs_numeric(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", ids[valid], ignore.case = TRUE))
  valid <- valid & is.finite(position) & position > 0 & position < 2^53
  data.frame(chromosome = chromosome, position = position, valid = valid)
}

minimum_cs_distance_assignment <- function(add_position, mix_position,
                                            add_id = seq_along(add_position),
                                            mix_id = seq_along(mix_position)) {
  k <- length(add_position)
  if (!k || length(mix_position) != k || length(add_id) != k || length(mix_id) != k ||
      any(!is.finite(c(add_position, mix_position)))) {
    stop("Matching requires equal, nonempty sets of finite positions and their IDs.")
  }
  a <- order(add_position, as.character(add_id), method = "radix")
  m <- order(mix_position, as.character(mix_id), method = "radix")
  distance <- abs(add_position[a] - mix_position[m])
  # A cost-neutral adjacent swap witnesses a nonunique minimum-total solution.
  # Keep this diagnostic separate from the deterministic secondary tie rule.
  nonunique <- FALSE
  if (k > 1L) {
    j <- seq_len(k - 1L)
    swapped <- abs(add_position[a[j]] - mix_position[m[j + 1L]]) +
      abs(add_position[a[j + 1L]] - mix_position[m[j]])
    nonunique <- any(swapped == distance[j] + distance[j + 1L])
  }
  list(add_index = a, mix_index = m, distance_bp = distance,
       minimum_total_distance_bp = sum(distance),
       minimum_total_has_tie = nonunique)
}

empty_matched_cs_pairs <- function() {
  data.frame(
    region_id = character(), gene = character(), tissue = character(),
    n_cs = integer(), pair_rank = integer(), region_weight = numeric(),
    chromosome = character(), genome_build = character(),
    additive_cs_id = character(), additive_cs_number = integer(),
    additive_cs_name = character(), additive_component = numeric(),
    additive_lead_snp = character(), additive_lead_predictor = character(),
    additive_lead_coding = character(), additive_lead_position = numeric(),
    additive_lead_pip = numeric(), additive_lead_alpha = numeric(),
    mixed_cs_id = character(), mixed_cs_number = integer(),
    mixed_cs_name = character(), mixed_component = numeric(),
    mixed_lead_snp = character(), mixed_lead_predictor = character(),
    mixed_lead_coding = character(), mixed_lead_position = numeric(),
    mixed_lead_pip = numeric(), mixed_lead_alpha = numeric(),
    distance_bp = numeric(), distance_kb = numeric(), zero_distance = logical(),
    same_lead_snp = logical(), scatter_eligible = logical(),
    region_total_distance_bp = numeric(), region_max_distance_bp = numeric(),
    minimum_total_has_tie = logical(), stringsAsFactors = FALSE
  )
}

match_equal_cs_regions <- function(res_summary, res_cs_summary) {
  regions <- as.data.frame(res_summary, stringsAsFactors = FALSE)
  cs <- as.data.frame(res_cs_summary, stringsAsFactors = FALSE)
  matched_cs_require(regions, c("gene", "tissue", "ncs_susie", "ncs_susie_mix"), "res_summary")
  cs_required <- c("gene", "tissue", "model_key", "cs_number", "lead_snp", "lead_coding", "lead_pip")
  if (!nrow(cs) && !length(names(cs))) {
    cs <- as.data.frame(setNames(rep(list(character()), length(cs_required)), cs_required))
  }
  matched_cs_require(cs, cs_required, "res_cs_summary")
  for (nm in c("gene", "tissue")) {
    regions[[nm]] <- as.character(regions[[nm]])
    cs[[nm]] <- as.character(cs[[nm]])
    if (anyNA(regions[[nm]]) || any(!nzchar(trimws(regions[[nm]]))) ||
        anyNA(cs[[nm]]) || any(!nzchar(trimws(cs[[nm]])))) {
      stop("Missing/empty ", nm, " identifiers; matching would be ambiguous.")
    }
  }
  if (anyDuplicated(regions[c("gene", "tissue")])) {
    stop("res_summary contains duplicate gene-tissue regions; resolve these before matching.")
  }
  regions <- regions[order(regions$gene, regions$tissue, method = "radix"), , drop = FALSE]
  rownames(regions) <- NULL
  regions$region_id <- matched_cs_key(regions$gene, regions$tissue)
  # Explicit model keys avoid inadvertently including weighted or permuted fits.
  cs <- cs[which(cs$model_key %in% c("susie_add", "susie_mix")), , drop = FALSE]
  cs$region_id <- matched_cs_key(cs$gene, cs$tissue)
  cs$cs_number <- matched_cs_numeric(cs$cs_number)
  for (nm in c("component", "lead_alpha")) {
    if (!nm %in% names(cs)) cs[[nm]] <- rep(NA_real_, nrow(cs))
  }
  for (nm in c("cs_name", "lead_predictor")) {
    if (!nm %in% names(cs)) cs[[nm]] <- rep(NA_character_, nrow(cs))
  }
  for (nm in c("lead_pip", "lead_alpha", "component")) cs[[nm]] <- matched_cs_numeric(cs[[nm]])
  cs$lead_coding <- tolower(as.character(cs$lead_coding))
  cs$lead_snp <- as.character(cs$lead_snp)
  # Build the identifier from the verified key and CS number, not an optional
  # preformatted identifier that might be duplicated or stale.
  cs$cs_id <- if (nrow(cs)) paste(cs$region_id, cs$model_key, cs$cs_number, sep = "|") else character()
  coordinates <- matched_cs_parse_snps(cs$lead_snp)
  cs$parsed_chromosome <- coordinates$chromosome
  cs$parsed_position <- coordinates$position
  cs$valid_coordinates <- coordinates$valid

  # Preserve orphan rows for diagnosis; they cannot supply a region count.
  orphan <- !cs$region_id %in% regions$region_id
  orphan_cs <- cs[orphan, , drop = FALSE]
  cs <- cs[!orphan, , drop = FALSE]
  by_region <- split(seq_len(nrow(cs)), factor(cs$region_id, levels = regions$region_id))
  n <- nrow(regions)
  audit <- regions
  audit$status <- rep("pending", n)
  audit$detail <- rep("", n)
  audit$available_additive_cs <- integer(n)
  audit$available_mixed_cs <- integer(n)
  audit$n_pairs <- integer(n)
  audit$n_scatter_points <- integer(n)
  audit$n_zero_distance <- integer(n)
  audit$minimum_total_distance_bp <- rep(NA_real_, n)
  audit$max_distance_bp <- rep(NA_real_, n)
  audit$minimum_total_has_tie <- rep(NA, n)
  rows <- vector("list", n)
  pair_columns <- names(empty_matched_cs_pairs())
  counts_add <- matched_cs_numeric(regions$ncs_susie)
  counts_mix <- matched_cs_numeric(regions$ncs_susie_mix)
  for (i in seq_len(n)) {
    part <- cs[by_region[[i]], , drop = FALSE]
    add <- part[part$model_key == "susie_add", , drop = FALSE]
    mix <- part[part$model_key == "susie_mix", , drop = FALSE]
    audit$available_additive_cs[i] <- nrow(add)
    audit$available_mixed_cs[i] <- nrow(mix)
    count <- c(counts_add[i], counts_mix[i])
    if (any(!is.finite(count)) || any(count < 0 | count != floor(count))) {
      audit$status[i] <- "invalid_cs_count"
      next
    }
    if (count[1] != count[2]) {
      audit$status[i] <- "unequal_cs_counts"
      next
    }
    k <- count[1]
    if (nrow(add) != k || nrow(mix) != k) {
      audit$status[i] <- "incomplete_or_inconsistent_cs_summary"
      audit$detail[i] <- "Observed CS row counts disagree with res_summary; no partial matching."
      next
    }
    if (!k) {
      audit$status[i] <- "equal_zero_cs"
      next
    }
    if (any(!is.finite(part$cs_number)) || any(part$cs_number != floor(part$cs_number)) ||
        !setequal(add$cs_number, seq_len(k)) || !setequal(mix$cs_number, seq_len(k))) {
      audit$status[i] <- "invalid_or_duplicate_cs_number"
      next
    }
    if (any(!part$valid_coordinates)) {
      audit$status[i] <- "invalid_lead_coordinates"
      audit$detail[i] <- paste(part$cs_id[!part$valid_coordinates], collapse = ";")
      next
    }
    if (length(unique(part$parsed_chromosome)) != 1L) {
      audit$status[i] <- "inconsistent_chromosomes"
      next
    }
    # A different results file within a region signals mismatched input runs.
    if ("result_file" %in% names(part) && "result_file" %in% names(regions)) {
      files <- as.character(c(regions$result_file[i], part$result_file))
      files <- files[!is.na(files) & nzchar(files)]
      if (length(unique(basename(gsub("\\\\", "/", files)))) > 1L) {
        audit$status[i] <- "inconsistent_result_files"
        next
      }
    }
    # Summary positional metadata can be reconstructed from workhorse/TSS
    # fields. Keep disagreements visible but always use the actual lead SNP.
    metadata_disagrees <- FALSE
    if ("lead_position" %in% names(part)) {
      old_position <- matched_cs_numeric(part$lead_position)
      metadata_disagrees <- any(is.finite(old_position) & old_position != part$parsed_position)
    }
    matching <- minimum_cs_distance_assignment(
      add$parsed_position, mix$parsed_position,
      paste(add$lead_snp, sprintf("%09d", as.integer(add$cs_number)), sep = "|"),
      paste(mix$lead_snp, sprintf("%09d", as.integer(mix$cs_number)), sep = "|")
    )
    add <- add[matching$add_index, , drop = FALSE]
    mix <- mix[matching$mix_index, , drop = FALSE]
    pair <- data.frame(
      region_id = regions$region_id[i], gene = regions$gene[i], tissue = regions$tissue[i],
      n_cs = k, pair_rank = seq_len(k), region_weight = 1 / k,
      chromosome = paste0("chr", add$parsed_chromosome), genome_build = "GRCh38",
      stringsAsFactors = FALSE
    )
    for (prefix in c("additive", "mixed")) {
      source <- if (prefix == "additive") add else mix
      for (nm in c("cs_id", "cs_number", "cs_name", "component", "lead_snp",
                   "lead_predictor", "lead_coding", "lead_pip", "lead_alpha")) {
        pair[[paste(prefix, nm, sep = "_")]] <- source[[nm]]
      }
      pair[[paste0(prefix, "_lead_position")]] <- source$parsed_position
    }
    pair$distance_bp <- matching$distance_bp
    pair$distance_kb <- pair$distance_bp / 1000
    pair$zero_distance <- pair$distance_bp == 0
    pair$same_lead_snp <- add$lead_snp == mix$lead_snp
    pair$scatter_eligible <- is.finite(mix$lead_pip) & mix$lead_pip >= 0 & mix$lead_pip <= 1
    pair$region_total_distance_bp <- matching$minimum_total_distance_bp
    pair$region_max_distance_bp <- max(pair$distance_bp)
    pair$minimum_total_has_tie <- matching$minimum_total_has_tie
    rows[[i]] <- pair[, pair_columns, drop = FALSE]
    audit$status[i] <- "matched"
    if (metadata_disagrees) audit$detail[i] <- "Saved lead_position differs; used coordinates parsed from lead_snp."
    audit$n_pairs[i] <- k
    audit$n_scatter_points[i] <- sum(pair$scatter_eligible)
    audit$n_zero_distance[i] <- sum(pair$zero_distance)
    audit$minimum_total_distance_bp[i] <- matching$minimum_total_distance_bp
    audit$max_distance_bp[i] <- max(pair$distance_bp)
    audit$minimum_total_has_tie[i] <- matching$minimum_total_has_tie
  }
  nonempty <- !vapply(rows, is.null, logical(1))
  pairs <- if (any(nonempty)) do.call(rbind, rows[nonempty]) else empty_matched_cs_pairs()
  # Copy scalar selection metadata to each pair for future tissue/cohort work.
  metadata <- intersect(c("result_file", "min_pv", "mean_count", "n_ind", "n_SNP",
                          "converged_add", "converged_mix", "log_lik_add", "log_lik_mix"), names(regions))
  for (nm in metadata) pairs[[nm]] <- regions[[nm]][match(pairs$region_id, regions$region_id)]
  rownames(pairs) <- NULL
  stopifnot(nrow(pairs) == sum(audit$n_pairs),
            !anyDuplicated(pairs[c("region_id", "additive_cs_id")]),
            !anyDuplicated(pairs[c("region_id", "mixed_cs_id")]))
  list(pairs = pairs, region_audit = audit, orphan_cs = orphan_cs)
}

matched_cs_distance_summary <- function(pairs, tissue = "All tissues") {
  distance <- pairs$distance_bp
  q <- if (length(distance)) as.numeric(quantile(distance, c(.25, .5, .75, .9, .95))) else rep(NA_real_, 5)
  data.frame(
    tissue = tissue, n_regions = length(unique(pairs$region_id)), n_pairs = nrow(pairs),
    n_zero_distance = sum(distance == 0), n_same_lead_snp = sum(pairs$same_lead_snp),
    n_positive_distance = sum(distance > 0), n_over_100kb = sum(distance > 100000),
    fraction_zero_distance = if (length(distance)) mean(distance == 0) else NA_real_,
    mean_distance_bp = if (length(distance)) mean(distance) else NA_real_,
    q25_distance_bp = q[1], median_distance_bp = q[2], q75_distance_bp = q[3],
    q90_distance_bp = q[4], q95_distance_bp = q[5],
    max_distance_bp = if (length(distance)) max(distance) else NA_real_,
    n_missing_or_invalid_mixed_pip = sum(!pairs$scatter_eligible)
  )
}

plot_matched_cs_distances <- function(pairs, output_file, tissue = NULL) {
  if (!is.null(tissue)) pairs <- pairs[which(pairs$tissue == tissue), , drop = FALSE]
  valid_pip <- is.finite(pairs$mixed_lead_pip) & pairs$mixed_lead_pip >= 0 & pairs$mixed_lead_pip <= 1
  shown <- pairs[valid_pip, , drop = FALSE]
  zero <- shown[shown$distance_bp == 0, , drop = FALSE]
  positive <- shown[shown$distance_bp > 0, , drop = FALSE]
  coding_names <- c("additive", "dominant", "recessive", "unknown")
  palette <- c(additive = "#68849c", dominant = "#087e83", recessive = "#d58046", unknown = "#777777")
  coding <- function(d) ifelse(d$mixed_lead_coding %in% names(palette), d$mixed_lead_coding, "unknown")
  draw_points <- function(d, x) {
    if (!nrow(d)) return(invisible(NULL))
    points(x, d$mixed_lead_pip, pch = ifelse(d$n_cs == 1L, 16, 4), cex = .53,
           col = grDevices::adjustcolor(palette[coding(d)], alpha.f = .38), lwd = .7)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  extension <- tolower(tools::file_ext(output_file))
  if (extension == "pdf") {
    grDevices::pdf(output_file, width = 10, height = 7.8, useDingbats = FALSE)
  } else if (extension == "png") {
    grDevices::png(output_file, width = 2200, height = 1716, res = 220,
                  type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  } else stop("Plot output must be .png or .pdf.")
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.35, 6))
  par(oma = c(3.7, 0, 6.7, .5), family = "sans", col.axis = "#444444", col.lab = "#333333")
  par(mar = c(4.4, 5.1, .7, .6))
  plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(-.025, 1.025), axes = FALSE,
       xlab = "", ylab = "SuSiE-mix lead PIP (lead coding)", yaxs = "i")
  axis(2, at = seq(0, 1, .25), las = 1)
  axis(1, at = .5, labels = "0 bp")
  abline(h = .5, col = "#d9d9d9", lty = 2)
  box(bty = "l", col = "#444444")
  # Deterministic horizontal spread is ONLY within the labelled zero strip;
  # it preserves each PIP, uses no RNG and never changes distance values.
  zero_x <- if (nrow(zero) <= 1L) rep(.5, nrow(zero)) else .15 + .7 * ((seq_len(nrow(zero)) * .61803398875) %% 1)
  draw_points(zero, zero_x)
  mtext(sprintf("(n = %s)", format(nrow(zero), big.mark = ",")), side = 1, line = 2.8, cex = .73)
  par(mar = c(4.4, .7, .7, 1.2))
  log_limits <- if (nrow(positive)) range(log10(positive$distance_kb)) else c(-3, 3)
  if (diff(log_limits) < .5) log_limits <- mean(log_limits) + c(-.35, .35)
  log_limits <- log_limits + c(-.09, .09) * diff(log_limits)
  plot(NA_real_, NA_real_, xlim = 10^log_limits, ylim = c(-.025, 1.025), log = "x", axes = FALSE,
       xlab = "Lead-SNP distance (kb, logarithmic scale)", ylab = "", yaxs = "i", xaxs = "i")
  ticks <- if (ceiling(log_limits[1]) <= floor(log_limits[2])) {
    10^seq(ceiling(log_limits[1]), floor(log_limits[2]))
  } else 10^mean(log_limits)
  tick_labels <- vapply(ticks, function(x) format(x, trim = TRUE, scientific = FALSE, big.mark = ","), character(1))
  axis(1, at = ticks, labels = tick_labels)
  abline(h = .5, col = "#d9d9d9", lty = 2)
  if (log_limits[1] < 2 && log_limits[2] > 2) abline(v = 100, col = "#d9d9d9", lty = 2)
  box(bty = "l", col = "#444444")
  draw_points(positive, positive$distance_kb)
  if (!nrow(positive)) text(10^mean(log_limits), .55, "No positive-distance pairs with a valid PIP", cex = .85)
  # Use a figure-wide overlay so legends do not depend on the log-axis range.
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(.065, .967, "Lead-SNP distance and support across matched credible sets", adj = 0, font = 2, cex = 1.1)
  cohort <- if (is.null(tissue)) "All tissues" else tissue
  text(.065, .929, sprintf("%s | %s regions | %s matched CS pairs | equal CS counts in both models",
                           cohort, format(length(unique(pairs$region_id)), big.mark = ","),
                           format(nrow(pairs), big.mark = ",")), adj = 0, cex = .83)
  counts <- table(factor(coding(shown), levels = coding_names))
  present <- counts > 0
  if (any(present)) legend(.065, .908, legend = paste0(c("Additive", "Dominant", "Recessive", "Unknown")[present],
      " (", format(as.integer(counts[present]), big.mark = ",", trim = TRUE), ")"),
      col = palette[present], pch = 16, horiz = TRUE, bty = "n", cex = .78, x.intersp = .6)
  legend(.065, .869, legend = c("1 CS per model", "2+ CSs per model"), pch = c(16, 4),
         col = "#555555", horiz = TRUE, bty = "n", cex = .76, x.intersp = .6)
  text(.065, .061, "One point per CS pair; one-to-one matching minimizes the total absolute GRCh38 distance.", adj = 0, cex = .75)
  text(.065, .036, sprintf("Zero distances use a separate strip. %s pairs lack a valid PIP; retained in the distance tables.",
                           sum(!valid_pip)), adj = 0, cex = .73, col = "#555555")
  invisible(list(n_pairs = nrow(pairs), n_zero = nrow(zero), n_positive = nrow(positive), n_missing_pip = sum(!valid_pip)))
}

read_matched_cs_input <- function(path, object_name) {
  if (!file.exists(path)) stop("Missing input: ", path, ". Run generate_summary_results.R first, or supply the summary paths.")
  if (tolower(tools::file_ext(path)) == "csv") {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) stop(path, " does not contain ", object_name, ".")
  as.data.frame(get(object_name, envir = env, inherits = FALSE))
}

run_matched_cs_lead_distances <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list("project-dir" = "/project2/mstephens/wdenault/susie_mix",
                   "summary-file" = NULL, "cs-summary-file" = NULL,
                   "output-dir" = NULL, "pairs-file" = NULL, "tissue" = NULL)
  if ("--help" %in% args) {
    cat("Match equal-count CS regions and plot lead distance against coding-specific mixed lead PIP.\n",
        "Usage: Rscript plot_matched_cs_lead_distances.R [--name=value ...]\n",
        "Options: ", paste(names(defaults), collapse = ", "), "\n",
        "Default inputs: PROJECT/res_summary.RData and PROJECT/res_cs_summary.RData.\n",
        "Default outputs: PROJECT/plot/matched_cs_lead_distances/.\n",
        "Use --pairs-file=... --tissue=Blood to replot saved pairs without rematching.\n", sep = "")
    return(invisible(NULL))
  }
  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) stop("Expected --name=value, got: ", arg)
    name <- sub("^--([^=]+)=.*$", "\\1", arg)
    if (!name %in% names(defaults)) stop("Unknown option: ", name)
    defaults[[name]] <- sub("^--[^=]+=", "", arg)
  }
  opt <- defaults
  project <- opt[["project-dir"]]
  tissue <- opt$tissue
  out <- opt[["output-dir"]]
  if (is.null(out)) {
    out <- if (is.null(opt[["pairs-file"]])) file.path(project, "plot", "matched_cs_lead_distances") else dirname(opt[["pairs-file"]])
    if (!is.null(tissue)) out <- file.path(out, "by_tissue", gsub("[^[:alnum:]_-]", "_", tissue))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (is.null(opt[["pairs-file"]])) {
    region_path <- if (is.null(opt[["summary-file"]])) file.path(project, "res_summary.RData") else opt[["summary-file"]]
    cs_path <- if (is.null(opt[["cs-summary-file"]])) file.path(project, "res_cs_summary.RData") else opt[["cs-summary-file"]]
    message("Reading region and CS summaries...")
    result <- match_equal_cs_regions(read_matched_cs_input(region_path, "res_summary"),
                                     read_matched_cs_input(cs_path, "res_cs_summary"))
    pairs <- result$pairs
    write.csv(pairs, file.path(out, "matched_cs_pairs.csv"), row.names = FALSE, na = "")
    write.csv(result$region_audit, file.path(out, "region_matching_audit.csv"), row.names = FALSE, na = "")
    write.csv(result$orphan_cs, file.path(out, "orphan_cs_rows.csv"), row.names = FALSE, na = "")
    print(table(result$region_audit$status))
    if (nrow(result$orphan_cs)) warning(nrow(result$orphan_cs), " CS rows have no region summary; see orphan_cs_rows.csv.")
    input_description <- paste("Region summary:", normalizePath(region_path, winslash = "/"),
                               "\nCS summary:", normalizePath(cs_path, winslash = "/"))
  } else {
    pairs <- read_matched_cs_input(opt[["pairs-file"]], "pairs")
    matched_cs_require(pairs, names(empty_matched_cs_pairs()), "matched pairs")
    if (any(!is.finite(pairs$distance_bp) | pairs$distance_bp < 0)) stop("Saved pair distances must be finite and nonnegative.")
    input_description <- paste("Previously matched pairs:", normalizePath(opt[["pairs-file"]], winslash = "/"))
  }
  if (!is.null(tissue) && !tissue %in% pairs$tissue) {
    stop("No matched pairs for tissue '", tissue, "'. Available: ", paste(sort(unique(pairs$tissue)), collapse = ", "))
  }
  selected <- if (is.null(tissue)) pairs else pairs[pairs$tissue == tissue, , drop = FALSE]
  write.csv(matched_cs_distance_summary(selected, if (is.null(tissue)) "All tissues" else tissue),
            file.path(out, "distance_summary.csv"), row.names = FALSE, na = "")
  # Save per-tissue descriptive summaries now, without making any independence
  # assumptions or testing tissues against an overlapping overall population.
  tissue_rows <- lapply(sort(unique(selected$tissue)), function(t) matched_cs_distance_summary(selected[selected$tissue == t, , drop = FALSE], t))
  per_tissue <- if (length(tissue_rows)) do.call(rbind, tissue_rows) else matched_cs_distance_summary(selected)[FALSE, ]
  write.csv(per_tissue, file.path(out, "distance_summary_by_tissue.csv"), row.names = FALSE, na = "")
  for (ext in c("png", "pdf")) plot_matched_cs_distances(pairs, file.path(out, paste0("matched_cs_lead_distances.", ext)), tissue)
  writeLines(c(
    "Equal-count SuSiE / unweighted SuSiE-mix lead distances", input_description,
    paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Plot tissue:", if (is.null(tissue)) "All tissues" else tissue),
    "Cohort: all input regions with equal positive CS counts; no association/read/convergence/coding filters.",
    "One row per matched CS pair. A k-CS region contributes k observations, including zero distances.",
    "Match objective: minimum total absolute physical distance; each CS is used once.",
    "Method: sort GRCh38 lead positions in each model, pair ranks. This is an exact optimum in one dimension.",
    "Tie rule: minimum sum of squared distances, then SNP identity / CS number. Does not maximize unchanged SNPs.",
    "minimum_total_has_tie flags a cost-neutral adjacent swap, i.e. a nonunique primary optimum.",
    "These positional pairs are not evidence of biological signal equivalence.",
    "Distance uses lead_snp coordinates, not distance to TSS. Incomplete or invalid regions are excluded in full and audited.",
    "same_lead_snp compares full SNP identifiers. zero_distance compares base-pair positions; these need not coincide.",
    "Y = mixed lead predictor PIP for its coding; neither CS alpha nor summed SNP PIP. Missing/invalid PIP affects scatter only.",
    "Every zero-distance point is spread horizontally within a separate strip; its PIP and saved distance are unchanged.",
    "Distance summaries include all matched pairs, including zeros and missing-PIP rows; quantiles use R type 7.",
    "Pairs are CS-weighted. region_weight = 1 / n_cs is supplied for a future equal-region-weighted comparison.",
    "Gene and tissue remain on each row. Pairs within a region and genes across tissues may be dependent.",
    "For a tissue-versus-background inferential comparison, define a nonoverlapping reference and account for clustering.",
    "matched_cs_pairs.csv contains all matched tissues, even when --tissue selects only the plot and summaries."
  ), file.path(out, "matching_notes.txt"))
  message("Saved ", nrow(selected), " pairs' distance/support plot and tables to ", normalizePath(out, winslash = "/"))
  invisible(list(pairs = pairs, output_dir = out))
}

if (sys.nframe() == 0L) run_matched_cs_lead_distances()
