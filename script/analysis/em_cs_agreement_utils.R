# Summary-only comparison of additive SuSiE and final EM weighted_fit_mix.
# A region is one gene/tissue. CS matching is descriptive, not evidence of LD
# or signal equivalence. Only base R is required.
em_agreement_require <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

em_agreement_key <- function(gene, tissue) {
  if (!length(gene)) return(character())
  paste0(nchar(gene), ":", gene, "|", nchar(tissue), ":", tissue)
}

em_agreement_coordinates <- function(ids) {
  valid <- !is.na(ids) & grepl("^chr([0-9]+|X|Y|M|MT)_[0-9]+_[^_]+_[^_]+_b38(_[^_]+)?$",
                             ids, ignore.case = TRUE)
  chr <- rep(NA_character_, length(ids))
  pos <- rep(NA_real_, length(ids))
  chr[valid] <- toupper(sub("^chr([^_]+)_.*$", "\\1", ids[valid], ignore.case = TRUE))
  chr[chr %in% "MT"] <- "M"
  pos[valid] <- as.numeric(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", ids[valid], ignore.case = TRUE))
  data.frame(chr = chr, pos = pos, valid = valid & is.finite(pos) & pos > 0 & pos < 2^53)
}

# Rectangular minimum-total-absolute-distance matching. Every element of the
# smaller set is used once; the larger set can skip elements. In one dimension
# an order-preserving optimum exists, so a small dynamic program suffices.
# Equal costs keep the earlier sorted choice. IDs make ties reproducible.
em_rectangular_assignment <- function(a, m, aid = seq_along(a), mid = seq_along(m)) {
  empty <- matrix(integer(), ncol = 2L, dimnames = list(NULL, c("add", "em")))
  if (!length(a) || !length(m)) return(empty)
  if (any(!is.finite(c(a, m)))) stop("Matching requires finite coordinates.")
  swapped <- length(a) > length(m)
  x <- if (swapped) m else a
  y <- if (swapped) a else m
  xo <- order(x, as.character(if (swapped) mid else aid), method = "radix")
  yo <- order(y, as.character(if (swapped) aid else mid), method = "radix")
  x <- x[xo]; y <- y[yo]
  n <- length(x); k <- length(y)
  cost <- matrix(Inf, n + 1L, k + 1L)
  take <- matrix(FALSE, n, k)
  cost[1L, ] <- 0
  for (i in seq_len(n)) for (j in seq_len(k)) {
    include <- cost[i, j] + abs(x[i] - y[j])
    skip <- cost[i + 1L, j]
    take[i, j] <- include < skip
    cost[i + 1L, j + 1L] <- min(include, skip)
  }
  selected <- matrix(NA_integer_, n, 2L)
  i <- n; j <- k
  while (i > 0L) {
    if (take[i, j]) {
      selected[i, ] <- c(xo[i], yo[j])
      i <- i - 1L
    }
    j <- j - 1L
  }
  if (swapped) selected <- selected[, 2:1, drop = FALSE]
  colnames(selected) <- c("add", "em")
  selected
}

# Reserve identical biological lead IDs first, ignoring predictor coding.
# Repeated leads are treated as a multiset: no CS is reused.
em_pair_cs_leads <- function(add, em) {
  pairs <- matrix(integer(), ncol = 2L, dimnames = list(NULL, c("add", "em")))
  for (snp in sort(intersect(add$lead_snp, em$lead_snp))) {
    a <- which(add$lead_snp == snp); m <- which(em$lead_snp == snp)
    a <- a[order(add$cs_number[a])]; m <- m[order(em$cs_number[m])]
    n <- min(length(a), length(m))
    pairs <- rbind(pairs, cbind(add = head(a, n), em = head(m, n)))
  }
  a <- setdiff(seq_len(nrow(add)), pairs[, "add"])
  m <- setdiff(seq_len(nrow(em)), pairs[, "em"])
  rest <- em_rectangular_assignment(add$position[a], em$position[m],
                                    add$cs_number[a], em$cs_number[m])
  if (nrow(rest)) pairs <- rbind(pairs, cbind(add = a[rest[, "add"]], em = m[rest[, "em"]]))
  pairs
}

summarize_em_cs_agreement <- function(res_summary, res_cs_summary,
    association_threshold = 1e-8, minimum_mean_reads = 100, minimum_tissue_regions = 100L) {
  if (length(association_threshold) != 1L || !is.finite(association_threshold) || association_threshold <= 0 ||
      length(minimum_mean_reads) != 1L || !is.finite(minimum_mean_reads) || minimum_mean_reads < 0 ||
      length(minimum_tissue_regions) != 1L || !is.finite(minimum_tissue_regions) ||
      minimum_tissue_regions < 1 || minimum_tissue_regions != floor(minimum_tissue_regions))
    stop("Invalid association, read-count, or minimum-tissue-size setting.")
  r <- as.data.frame(res_summary, stringsAsFactors = FALSE)
  cs <- as.data.frame(res_cs_summary, stringsAsFactors = FALSE)
  em_agreement_require(r, c("gene", "tissue", "min_pv", "mean_count", "ncs_susie",
    "ncs_weighted_fit_mix", "has_weighted_fit_mix", "em_iteration"), "EM region summary")
  em_agreement_require(cs, c("gene", "tissue", "model_key", "cs_number", "lead_snp",
    "lead_coding", "lead_pip", "em_iteration"), "EM CS summary")
  # Avoid copying large CS-membership strings and unrelated summary metrics
  # into the audit tables. Counts here have their own, explicit meanings.
  r <- r[c("gene", "tissue", "min_pv", "mean_count", "ncs_susie", "ncs_weighted_fit_mix",
    "has_weighted_fit_mix", "em_iteration", intersect(c("result_file", "n_ind",
    "converged_add", "converged_weighted_fit_mix"), names(r)))]
  cs <- cs[c("gene", "tissue", "model_key", "cs_number", "lead_snp", "lead_coding", "lead_pip",
    "em_iteration", intersect(c("component", "cs_size_snps"), names(cs)))]
  iteration <- unique(r$em_iteration)
  if (length(iteration) != 1L || !is.finite(iteration) ||
      anyNA(cs$em_iteration) || any(cs$em_iteration != iteration)) stop("Expected one EM iteration in both summaries.")
  for (nm in c("gene", "tissue")) {
    r[[nm]] <- as.character(r[[nm]]); cs[[nm]] <- as.character(cs[[nm]])
    if (anyNA(r[[nm]]) || any(!nzchar(r[[nm]])) || anyNA(cs[[nm]]) || any(!nzchar(cs[[nm]])))
      stop("Missing gene/tissue identifiers.")
  }
  if (anyDuplicated(r[c("gene", "tissue")])) stop("Duplicate gene/tissue regions.")
  r <- r[order(r$tissue, r$gene), , drop = FALSE]
  rownames(r) <- NULL
  r$region_id <- em_agreement_key(r$gene, r$tissue)
  r$n_add <- as.numeric(r$ncs_susie); r$n_em <- as.numeric(r$ncs_weighted_fit_mix)
  valid_count <- is.finite(r$n_add) & is.finite(r$n_em) & r$n_add >= 0 & r$n_em >= 0 &
    r$n_add == floor(r$n_add) & r$n_em == floor(r$n_em)
  r$selection_status <- "included"
  r$selection_status[which(valid_count & r$n_add + r$n_em == 0)] <- "both_zero_cs"
  r$selection_status[!valid_count] <- "invalid_cs_count"
  r$selection_status[!is.finite(r$mean_count) | r$mean_count < minimum_mean_reads] <- "read_count_filter"
  r$selection_status[!is.finite(r$min_pv) | r$min_pv >= association_threshold] <- "association_filter"
  r$selection_status[!r$has_weighted_fit_mix %in% TRUE] <- "missing_em_fit"
  r$count_eligible <- r$selection_status == "included"
  r$count_group <- ifelse(!valid_count, NA_character_, ifelse(r$n_add == r$n_em,
    "Equal CS counts", ifelse(r$n_em > r$n_add, "More EM CSs", "More SuSiE CSs")))
  r$lead_status <- r$matching_status <- "not_selected"
  r$n_shared_leads <- r$n_pairs <- r$n_extra_add <- r$n_extra_em <- NA_integer_
  r$retention <- NA_character_
  cs <- cs[which(cs$model_key %in% c("susie_add", "weighted_fit_mix")), , drop = FALSE]
  cs$region_id <- em_agreement_key(cs$gene, cs$tissue)
  cs$lead_snp <- as.character(cs$lead_snp)
  cs$cs_id <- paste(cs$region_id, cs$model_key, cs$cs_number, sep = "|")
  coordinates <- em_agreement_coordinates(cs$lead_snp)
  cs$position <- coordinates$pos; cs$chromosome <- coordinates$chr
  cs$valid_coordinates <- coordinates$valid
  orphan <- cs[!cs$region_id %in% r$region_id, , drop = FALSE]
  parts <- split(seq_len(nrow(cs)), factor(cs$region_id, levels = r$region_id))
  pair_rows <- extra_rows <- list()
  for (i in which(r$count_eligible)) {
    part <- cs[parts[[i]], , drop = FALSE]
    add <- part[part$model_key == "susie_add", , drop = FALSE]
    em <- part[part$model_key == "weighted_fit_mix", , drop = FALSE]
    invalid <- if (nrow(add) != r$n_add[i] || nrow(em) != r$n_em[i]) "incomplete_cs_summary" else
      if (anyNA(part$cs_number) || anyDuplicated(part$cs_id) ||
          !setequal(add$cs_number, seq_len(nrow(add))) || !setequal(em$cs_number, seq_len(nrow(em))))
        "invalid_cs_numbers" else if (anyNA(part$lead_snp) || any(!nzchar(part$lead_snp))) "missing_lead_id" else NULL
    if (!is.null(invalid)) {
      r$lead_status[i] <- r$matching_status[i] <- invalid
      next
    }
    r$lead_status[i] <- "ok"
    common <- intersect(add$lead_snp, em$lead_snp)
    n_shared <- sum(vapply(common, function(s) min(sum(add$lead_snp == s), sum(em$lead_snp == s)), integer(1)))
    r$n_shared_leads[i] <- n_shared
    smaller <- min(nrow(add), nrow(em))
    r$retention[i] <- if (!smaller) "Other model has no CS" else if (n_shared == smaller)
      "All leads retained" else if (n_shared > 0) "Some leads retained" else "No leads retained"
    # Identity-based retention remains usable if distances cannot be computed.
    if (any(!part$valid_coordinates) || length(unique(part$chromosome)) != 1L) {
      r$matching_status[i] <- "invalid_or_inconsistent_coordinates"
      next
    }
    r$matching_status[i] <- "matched"
    matched <- em_pair_cs_leads(add, em)
    r$n_pairs[i] <- nrow(matched)
    r$n_extra_add[i] <- nrow(add) - nrow(matched)
    r$n_extra_em[i] <- nrow(em) - nrow(matched)
    if (nrow(matched)) {
      a <- add[matched[, "add"], , drop = FALSE]; m <- em[matched[, "em"], , drop = FALSE]
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(region_id = r$region_id[i], gene = r$gene[i],
        tissue = r$tissue[i], n_add = nrow(add), n_em = nrow(em),
        additive_cs_id = a$cs_id, em_cs_id = m$cs_id, additive_lead_snp = a$lead_snp,
        em_lead_snp = m$lead_snp, em_lead_coding = m$lead_coding, em_lead_pip = m$lead_pip,
        same_lead_snp = a$lead_snp == m$lead_snp, distance_bp = abs(a$position - m$position))
    }
    for (model in c("SuSiE", "EM SuSiE-mix")) {
      own <- if (model == "SuSiE") add else em
      other <- if (model == "SuSiE") em else add
      used <- matched[, if (model == "SuSiE") "add" else "em"]
      extra <- own[setdiff(seq_len(nrow(own)), used), , drop = FALSE]
      if (!nrow(extra)) next
      nearest <- vapply(extra$position, function(pos) if (nrow(other)) min(abs(pos - other$position)) else NA_real_, numeric(1))
      extra_rows[[length(extra_rows) + 1L]] <- data.frame(region_id = r$region_id[i], gene = r$gene[i],
        tissue = r$tissue[i], model = model, n_add = nrow(add), n_em = nrow(em),
        cs_id = extra$cs_id, lead_snp = extra$lead_snp, lead_coding = extra$lead_coding,
        lead_pip = extra$lead_pip, nearest_other_distance_bp = nearest,
        other_model_has_cs = nrow(other) > 0)
    }
  }
  pairs <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame(
    region_id = character(), gene = character(), tissue = character(), n_add = integer(), n_em = integer(),
    additive_cs_id = character(), em_cs_id = character(), additive_lead_snp = character(),
    em_lead_snp = character(), em_lead_coding = character(), em_lead_pip = numeric(),
    same_lead_snp = logical(), distance_bp = numeric())
  extras <- if (length(extra_rows)) do.call(rbind, extra_rows) else data.frame(
    region_id = character(), gene = character(), tissue = character(), model = character(),
    n_add = integer(), n_em = integer(), cs_id = character(), lead_snp = character(),
    lead_coding = character(), lead_pip = numeric(), nearest_other_distance_bp = numeric(), other_model_has_cs = logical())
  selected_regions <- r[r$count_eligible, , drop = FALSE]
  tissue_rows <- lapply(sort(unique(r$tissue)), function(tissue) {
    d <- selected_regions[selected_regions$tissue == tissue, , drop = FALSE]
    same <- sum(d$n_add == d$n_em)
    data.frame(tissue = tissue, n_regions = nrow(d), same_cs_count = same,
      different_cs_count = nrow(d) - same, more_em_cs = sum(d$n_em > d$n_add),
      more_additive_cs = sum(d$n_add > d$n_em),
      agreement_rate = if (nrow(d)) same / nrow(d) else NA_real_,
      n_both_zero_excluded = sum(r$tissue == tissue & r$selection_status == "both_zero_cs"),
      n_valid_lead_regions = sum(d$lead_status == "ok"), n_distance_regions = sum(d$matching_status == "matched"))
  })
  tissues <- do.call(rbind, tissue_rows)
  tissues$ranking_eligible <- tissues$n_regions >= minimum_tissue_regions
  ranked <- tissues[which(tissues$ranking_eligible), , drop = FALSE]
  low <- head(ranked[order(ranked$agreement_rate, -ranked$n_regions, ranked$tissue), , drop = FALSE],
              min(2L, max(1L, nrow(ranked) %/% 2L)))
  high <- ranked[!ranked$tissue %in% low$tissue, , drop = FALSE]
  high <- head(high[order(-high$agreement_rate, -high$n_regions, high$tissue), , drop = FALSE], 2L)
  low$selection <- rep("Most disagreement", nrow(low)); high$selection <- rep("Most agreement", nrow(high))
  featured <- rbind(low, high)
  if (nrow(featured) < 4L) warning("Only ", nrow(featured), " tissues meet minimum_tissue_regions = ",
                                 minimum_tissue_regions, "; fewer than four tissue tables will be shown.", call. = FALSE)
  count_tables <- lapply(seq_len(nrow(featured)), function(i) {
    d <- selected_regions[selected_regions$tissue == featured$tissue[i], , drop = FALSE]
    levels <- sort(unique(c(0L, d$n_add, d$n_em)))
    tab <- as.data.frame(table(SuSiE = factor(d$n_add, levels), EM = factor(d$n_em, levels)), stringsAsFactors = FALSE)
    data.frame(tissue = featured$tissue[i], additive_cs_count = as.integer(as.character(tab$SuSiE)),
      em_cs_count = as.integer(as.character(tab$EM)), n_regions = tab$Freq)
  })
  tables <- if (length(count_tables)) do.call(rbind, count_tables) else data.frame(
    tissue = character(), additive_cs_count = integer(), em_cs_count = integer(), n_regions = integer())
  group_levels <- c("Equal CS counts", "More EM CSs", "More SuSiE CSs")
  retention_levels <- c("All leads retained", "Some leads retained", "No leads retained", "Other model has no CS")
  retained <- r[r$count_eligible & r$lead_status == "ok", , drop = FALSE]
  retention <- as.data.frame(table(count_group = factor(retained$count_group, group_levels),
    retention = factor(retained$retention, retention_levels)), stringsAsFactors = FALSE)
  names(retention)[3] <- "n_regions"
  retention$denominator <- ave(retention$n_regions, retention$count_group, FUN = sum)
  retention$proportion <- ifelse(retention$denominator > 0, retention$n_regions / retention$denominator, NA_real_)
  distances <- list("Paired, changed lead" = pairs$distance_bp[!pairs$same_lead_snp],
    "Extra EM to nearest SuSiE" = extras$nearest_other_distance_bp[extras$model == "EM SuSiE-mix"],
    "Extra SuSiE to nearest EM" = extras$nearest_other_distance_bp[extras$model == "SuSiE"])
  distance_summary <- do.call(rbind, lapply(names(distances), function(group) {
    raw <- distances[[group]]; d <- raw[is.finite(raw)]
    data.frame(group = group, n_total = length(raw), n_distance = length(d), n_without_counterpart = sum(!is.finite(raw)),
      n_zero_distance = sum(d == 0), median_kb = if (length(d)) median(d) / 1000 else NA_real_,
      p90_kb = if (length(d)) unname(quantile(d, .9)) / 1000 else NA_real_,
      n_over_10kb = sum(d > 1e4), n_over_100kb = sum(d > 1e5))
  }))
  extra_em <- extras[extras$model == "EM SuSiE-mix", , drop = FALSE]
  extra_em$coding_group <- ifelse(extra_em$lead_coding %in% c("additive", "recessive", "dominant"),
                                 as.character(extra_em$lead_coding), "unknown")
  codings <- c("additive", "recessive", "dominant", if (any(extra_em$coding_group == "unknown")) "unknown")
  coding_summary <- do.call(rbind, lapply(codings, function(coding) {
    d <- extra_em[extra_em$coding_group == coding, , drop = FALSE]
    pip <- d$lead_pip[is.finite(d$lead_pip) & d$lead_pip >= 0 & d$lead_pip <= 1]
    data.frame(coding = coding, n_extra_cs = nrow(d), n_valid_pip = length(pip),
      n_without_additive_cs = sum(!d$other_model_has_cs),
      median_pip = if (length(pip)) median(pip) else NA_real_, n_pip_at_least_half = sum(pip >= .5))
  }))
  stopifnot(!anyDuplicated(pairs$additive_cs_id), !anyDuplicated(pairs$em_cs_id),
            !any(extras$cs_id %in% c(pairs$additive_cs_id, pairs$em_cs_id)),
            nrow(pairs) == sum(r$n_pairs, na.rm = TRUE))
  list(iteration = iteration, region_audit = r, tissue_ranking = tissues, selected_tissues = featured,
    tissue_cs_count_tables = tables, matched_cs_pairs = pairs, unpaired_cs = extras,
    lead_retention_summary = retention, lead_distance_summary = distance_summary,
    extra_em_coding_summary = coding_summary, extra_em = extra_em, distances = distances,
    orphan_cs_rows = orphan, settings = list(association_threshold = association_threshold,
      minimum_mean_reads = minimum_mean_reads, minimum_tissue_regions = minimum_tissue_regions))
}

plot_em_cs_agreement <- function(result, output_file = NULL, preview_label = NULL) {
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (tolower(tools::file_ext(output_file)) == "pdf") {
      grDevices::pdf(output_file, width = 16, height = 12, useDingbats = FALSE)
    } else grDevices::png(output_file, width = 3200, height = 2400, res = 200,
                          type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE, after = FALSE)
  par(mfrow = c(2, 2), oma = c(2.8, .5, 4, .5), family = "sans", cex = .9,
      col.axis = "#444444", col.lab = "#333333")
  # A: four literal contingency tables; no binning or rounding of CS counts.
  par(mar = c(1, 1.5, 3, 1))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  title("A   CS-count agreement: four tissue extremes", adj = 0, cex.main = 1.1)
  text(.5, .99, "Rows: SuSiE CSs   |   Columns: EM CSs   |   Cells: region counts", cex = .78)
  if (!nrow(result$selected_tissues)) text(.5, .5, "No tissues meet the ranking sample-size threshold.")
  for (i in seq_len(nrow(result$selected_tissues))) {
    f <- result$selected_tissues[i, ]; tab <- result$tissue_cs_count_tables
    tab <- tab[tab$tissue == f$tissue, ]
    levels <- sort(unique(tab$additive_cs_count)); k <- length(levels)
    left <- if (i %% 2L) .055 else .56
    top <- if (i <= 2L) .79 else .33
    w <- .39 / k; h <- .29 / k
    heading_col <- if (f$selection == "Most disagreement") "#b75c32" else "#287d72"
    text(left + .195, top + .125, f$tissue, font = 2, cex = .9, col = heading_col)
    text(left + .195, top + .086, sprintf("%s | n = %s", f$selection, format(f$n_regions, big.mark = ",")), cex = .7)
    text(left + .195, top + .05, sprintf("Same: %d (%.1f%%)  |  Different: %d", f$same_cs_count,
      100 * f$agreement_rate, f$different_cs_count), cex = .7)
    for (a in seq_len(k)) for (m in seq_len(k)) {
      n <- tab$n_regions[tab$additive_cs_count == levels[a] & tab$em_cs_count == levels[m]]
      x <- left + (m - 1L) * w; y <- top - a * h
      base <- if (a == m) "#287d72" else "#c16a40"
      fill <- grDevices::adjustcolor(base, alpha.f = if (n) .10 + .58 * sqrt(n / max(tab$n_regions)) else .025)
      rect(x, y, x + w, y + h, col = fill, border = "#e0e4e7", lwd = .5)
      if (levels[a] == 0L && levels[m] == 0L) {
        text(x + w / 2, y + h / 2, "-", cex = min(.8, 7 / k), col = "#888888")
      } else text(x + w / 2, y + h / 2, n, cex = min(.8, 7 / k))
    }
    text(left - .015, top - (seq_len(k) - .5) * h, levels, adj = 1, cex = min(.8, 7 / k))
    text(left + (seq_len(k) - .5) * w, top + .014, levels, cex = min(.8, 7 / k))
  }
  # B: separate one-to-one pair distances from nearest-neighbour references.
  par(mar = c(5.4, 4.5, 3, 1.5))
  colors <- c("#346b94", "#b75c32", "#7b68a5")
  vals <- unlist(result$distances); vals <- vals[is.finite(vals)]
  xmax <- max(6, if (length(vals)) log10(1 + max(vals)) else 0)
  plot(NA, xlim = c(0, xmax), ylim = c(0, 1.03), axes = FALSE,
       xlab = "", ylab = "Cumulative fraction of CS leads")
  axis(1, at = log10(1 + c(0, 10^(0:6))), labels = c("0", "1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb"), cex.axis = .8)
  axis(2, at = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"), las = 1)
  abline(h = seq(0, 1, .25), col = "#eeeeee"); abline(v = log10(1 + c(1e4, 1e5)), lty = 3, col = "#cccccc")
  title("B   How far do leads differ? All tissues", adj = 0, cex.main = 1.1)
  labels <- character()
  for (j in seq_along(result$distances)) {
    d <- result$distances[[j]]; d <- sort(d[is.finite(d)])
    if (length(d)) {
      x <- sort(unique(d)); y <- findInterval(x, d) / length(d)
      lines(c(0, log10(1 + x), xmax), c(0, y, 1), type = "s", lwd = 2, col = colors[j], lty = j)
    }
    labels <- c(labels, sprintf("%s (n = %s)", names(result$distances)[j], format(length(d), big.mark = ",")))
  }
  legend("bottomright", labels, col = colors, lty = 1:3, lwd = 2, bty = "o", box.col = NA, bg = "white", cex = .77)
  if (!length(vals)) text(xmax / 2, .55, "No comparable non-identical or extra leads")
  mtext("Lead-SNP distance (logarithmic spacing; zero retained)", side = 1, line = 2.5, cex = .9)
  mtext(sprintf("Excluded here: %s identical pairs; %s extra CSs without a counterpart model CS",
    format(sum(result$matched_cs_pairs$same_lead_snp), big.mark = ","),
    format(sum(!result$unpaired_cs$other_model_has_cs), big.mark = ",")), side = 1, line = 4.2, cex = .68)
  # C: regions, not individual CSs; each count group has its own denominator.
  par(mar = c(5.5, 9, 3, 1.5))
  groups <- c("Equal CS counts", "More EM CSs", "More SuSiE CSs")
  cats <- c("All leads retained", "Some leads retained", "No leads retained", "Other model has no CS")
  fills <- c("#287d72", "#85b9a5", "#c16a40", "#c6cbd0")
  rt <- result$lead_retention_summary
  plot(NA, xlim = c(0, 1), ylim = c(.5, 4.2), axes = FALSE,
       xlab = "Fraction of regions within each CS-count group", ylab = "")
  axis(1, at = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"))
  title("C   Lead retention across CS-count groups", adj = 0, cex.main = 1.1)
  for (i in seq_along(groups)) {
    d <- rt[rt$count_group == groups[i], ]; d <- d[match(cats, d$retention), ]
    y <- 4 - i; n <- sum(d$n_regions)
    axis(2, at = y, labels = sprintf("%s\n(n = %s)", groups[i], format(n, big.mark = ",")), las = 1, tick = FALSE, cex.axis = .85)
    start <- 0
    for (j in seq_along(cats)) {
      p <- if (n) d$n_regions[j] / n else 0
      rect(start, y - .30, start + p, y + .30, col = fills[j], border = "white")
      if (p >= .09) text(start + p / 2, y, sprintf("%.0f%%", 100 * p), cex = .85,
                          col = if (j %in% c(1L, 3L)) "white" else "#222222")
      start <- start + p
    }
    if (!n) text(.5, y, "No eligible regions", col = "#777777")
  }
  legend("top", legend = cats, fill = fills, ncol = 2, bty = "n", cex = .72)
  # D: coding-specific PIP of the extra EM CS leads, including EM-only regions.
  par(mar = c(5.5, 4.5, 3, 1.5))
  codings <- result$extra_em_coding_summary$coding
  plot(NA, xlim = c(.5, length(codings) + .5), ylim = c(0, 1.04), xaxt = "n", las = 1,
       xlab = "", ylab = "EM lead PIP (coding-specific)")
  title("D   Coding and support of extra EM CSs", adj = 0, cex.main = 1.1)
  abline(h = .5, lty = 3, col = "#bbbbbb")
  coding_colors <- c(additive = "#346b94", recessive = "#9a6bb0", dominant = "#c18b35", unknown = "#999999")
  for (i in seq_along(codings)) {
    d <- result$extra_em[result$extra_em$coding_group == codings[i], ]
    p <- d$lead_pip[is.finite(d$lead_pip) & d$lead_pip >= 0 & d$lead_pip <= 1]
    if (length(p)) {
      boxplot(p, at = i, add = TRUE, axes = FALSE, boxwex = .45, outline = FALSE,
              col = grDevices::adjustcolor(coding_colors[codings[i]], .22), border = coding_colors[codings[i]])
      # Deterministic jitter, independent of the session random-number state.
      jitter <- ((seq_along(p) * .61803398875) %% 1 - .5) * .30
      points(i + jitter, p, pch = 16, cex = .35, col = grDevices::adjustcolor(coding_colors[codings[i]], .35))
    }
    axis(1, at = i, labels = sprintf("%s\nPIP n = %s / %s", tools::toTitleCase(codings[i]),
      format(length(p), big.mark = ","), format(nrow(d), big.mark = ",")), tick = FALSE, cex.axis = .85)
  }
  if (!nrow(result$extra_em)) text(mean(par("usr")[1:2]), .6, "No extra EM CSs")
  mtext("Extra = unpaired after identical-lead and distance matching; includes EM-only regions.", side = 1, line = 4, cex = .65)
  label <- paste0(if (!is.null(preview_label)) paste0(preview_label, " | ") else "",
                  "SuSiE versus EM SuSiE-mix | iteration ", result$iteration)
  mtext(label, outer = TRUE, side = 3, line = 2, cex = 1.35, font = 2)
  s <- result$settings; audit <- result$region_audit
  mtext(sprintf("P < %s | mean reads >= %s | at least one CS in either model | %s regions across %s tissues",
    format(s$association_threshold, scientific = TRUE), s$minimum_mean_reads,
    format(sum(audit$count_eligible), big.mark = ","), sum(result$tissue_ranking$n_regions > 0)),
    outer = TRUE, side = 3, line = .55, cex = .9, col = "#555555")
  mtext(sprintf("Tissue ranking: equal-count rate, minimum n = %d. CS leads compare exact SNP IDs, ignoring coding.",
    s$minimum_tissue_regions), outer = TRUE, side = 1, line = .7, cex = .78, col = "#555555")
  mtext(sprintf("Positional matches do not establish shared signals or LD. Audited exclusions: %d lead comparisons; %d distance comparisons.",
    sum(audit$count_eligible & audit$lead_status != "ok"), sum(audit$count_eligible & audit$matching_status != "matched")),
    outer = TRUE, side = 1, line = 1.9, cex = .78, col = "#555555")
  invisible(result)
}

write_em_cs_agreement <- function(result, output_dir, show_plot = interactive(), preview_label = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c("region_audit", "tissue_ranking", "selected_tissues", "tissue_cs_count_tables",
    "matched_cs_pairs", "unpaired_cs", "lead_retention_summary", "lead_distance_summary",
    "extra_em_coding_summary", "orphan_cs_rows")
  for (name in tables) write.csv(result[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  write.csv(data.frame(em_iteration = result$iteration, as.data.frame(result$settings)),
            file.path(output_dir, "analysis_settings.csv"), row.names = FALSE)
  for (ext in c("png", "pdf")) plot_em_cs_agreement(result,
    file.path(output_dir, paste0("em_cs_agreement_4panel.", ext)), preview_label)
  message("Selected tissue CS-count tables (rows = SuSiE; columns = EM):")
  for (tissue in result$selected_tissues$tissue) {
    cat("\n", tissue, "\n", sep = "")
    d <- result$tissue_cs_count_tables[result$tissue_cs_count_tables$tissue == tissue, ]
    print(xtabs(n_regions ~ additive_cs_count + em_cs_count, d))
  }
  message("Four-panel figure and CSV tables: ", normalizePath(output_dir, winslash = "/"))
  if (show_plot) plot_em_cs_agreement(result, preview_label = preview_label)
  invisible(result)
}
