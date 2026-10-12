# GWAS follow-up of a frozen, completed EM screening run.
cg_need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
}
cg_lse <- function(x) {
  m <- max(x)
  if (m == -Inf) return(-Inf)
  m + log(sum(exp(x - m)))
}
cg_key <- function(chromosome, position, a, b) {
  paste(sub("^chr", "", chromosome), format(position, scientific = FALSE, trim = TRUE),
        pmin(toupper(a), toupper(b)), pmax(toupper(a), toupper(b)), sep = ":")
}
cg_write <- function(rows, file) {
  x <- if (length(rows)) data.table::rbindlist(rows, fill = TRUE) else data.table::data.table(status = character())
  data.table::fwrite(x, file, na = "NA")
  invisible(x)
}

prepare_em_gwas_coloc <- function(project_dir, iteration = "latest", screen_dir = NULL,
                                  output_dir = NULL, padding_bp = 100000) {
  cg_need(c("jsonlite", "digest", "data.table"))
  context <- em_analysis_context(project_dir, iteration)
  if (is.null(screen_dir)) {
    dirs <- list.dirs(file.path(context$iteration_dir, "coloc_screen"), recursive = FALSE, full.names = TRUE)
    records <- lapply(dirs, function(d) tryCatch(readRDS(file.path(d, "screen_metadata.rds")), error = function(e) NULL))
    good <- vapply(records, function(m) !is.null(m$completed_at) && identical(as.integer(m$iteration), context$iteration), logical(1))
    if (!any(good)) stop("Run screen_em_coloc_candidates.R for the latest completed iteration first.")
    dirs <- dirs[good]; records <- records[good]
    screen_dir <- dirs[order(vapply(records, `[[`, "", "completed_at"), dirs, decreasing = TRUE)[1]]
  }
  screen_dir <- normalizePath(screen_dir, mustWork = TRUE)
  screen_meta <- readRDS(file.path(screen_dir, "screen_metadata.rds"))
  if (is.null(screen_meta$completed_at) || screen_meta$iteration != context$iteration)
    stop("Screen is incomplete or does not belong to the selected completed EM iteration.")
  manifest <- read.csv(file.path(screen_dir, "bf_manifest.csv"), stringsAsFactors = FALSE)
  manifest <- manifest[grepl("^exported", manifest$bf_status) & !is.na(manifest$bundle), , drop = FALSE]
  if (!nrow(manifest)) stop("No exported screening bundles.")
  if (is.null(output_dir)) output_dir <- tempfile(paste0("gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                                                tmpdir = file.path(context$iteration_dir, "coloc_gwas"))
  if (dir.exists(output_dir)) stop("Use a new output directory; existing jobs are frozen.")
  dir.create(file.path(output_dir, "bundles"), recursive = TRUE)
  dir.create(file.path(output_dir, "gwas"))
  regions <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    bundle <- readRDS(file.path(screen_dir, manifest$bundle[i]))
    if (bundle$schema_version != 1L || bundle$metadata$iteration != context$iteration) stop("Unexpected BF bundle version/iteration")
    v <- rbind(bundle$em$variants, bundle$vanilla$variants)
    chr <- unique(v$chromosome)
    if (length(chr) != 1L || !all(v$build == "GRCh38")) stop("Invalid bundle coordinates")
    id <- sprintf("case_%06d", i)
    path <- file.path("bundles", paste0(id, ".rds"))
    saveRDS(bundle, file.path(output_dir, path))
    manifest$case_id[i] <- id
    manifest$bundle[i] <- path
    manifest$md5[i] <- unname(tools::md5sum(file.path(output_dir, path)))
    regions[[i]] <- data.frame(case_id = id, chromosome = chr,
      start = max(1, min(v$position) - padding_bp), end = max(v$position) + padding_bp)
  }
  write.csv(manifest, file.path(output_dir, "cases.csv"), row.names = FALSE)
  write.csv(do.call(rbind, regions), file.path(output_dir, "regions.csv"), row.names = FALSE)
  file.copy(file.path(project_dir, "script/analysis/coloc_gwas_sources.json"), file.path(output_dir, "sources.json"))
  saveRDS(list(schema_version = 1L, iteration = context$iteration, screen_dir = screen_dir,
    screen_metadata = screen_meta, padding_bp = padding_bp,
    prepared_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(output_dir, "job_metadata.rds"))
  message("Frozen GWAS job: ", normalizePath(output_dir))
  normalizePath(output_dir)
}

cg_clean_gwas <- function(g, min_info = .8) {
  g <- as.data.frame(g)
  g$key <- cg_key(g$chromosome, g$position, g$effect_allele, g$other_allele)
  # Omit all duplicate allele keys; never select the smallest P value among duplicates.
  duplicate <- duplicated(g$key) | duplicated(g$key, fromLast = TRUE)
  pal <- paste0(pmin(g$effect_allele, g$other_allele), pmax(g$effect_allele, g$other_allele)) %in% c("AT", "CG")
  keep <- !duplicate & !pal & (is.na(g$info) | g$info >= min_info)
  attr(g, "qc") <- c(input = nrow(g), duplicate = sum(duplicate), palindromic = sum(pal), low_info = sum(!is.na(g$info) & g$info < min_info))
  list(data = g[keep, , drop = FALSE], counts = attr(g, "qc"))
}

cg_abf <- function(g, source, prior_sd = NULL) {
  if (!is.null(source$columns$log_bf)) {
    bf <- g$log_bf
    origin <- "published_MR_MEGA_logBF_single_signal_sensitivity"
    sd_used <- NA_real_
  } else {
    if (is.null(prior_sd)) {
      if (source$trait_type == "cc") prior_sd <- .2 else {
        sd_y <- source$sd_y
        if (is.null(sd_y)) {
          # coloc's variance relationship, using available per-SNP N and raw-scale SE.
          valid <- is.finite(g$n) & g$n > 0 & is.finite(g$eaf) & g$eaf > .01 & g$eaf < .99 &
            is.finite(g$se) & g$se > 0
          if (sum(valid) < 100L) stop("ABF needs a documented GWAS effect prior SD or sd_y; insufficient data to estimate sdY.")
          inv_v <- 1 / g$se[valid]^2
          x <- 2 * g$n[valid] * g$eaf[valid] * (1 - g$eaf[valid])
          sd_y <- sqrt(sum(inv_v * x) / sum(inv_v^2))
        }
        prior_sd <- .15 * sd_y
      }
    }
    if (length(prior_sd) != 1L || !is.finite(prior_sd) || prior_sd <= 0) stop("Invalid GWAS effect prior SD")
    r <- prior_sd^2 / (prior_sd^2 + g$se^2)
    bf <- .5 * (log1p(-r) + r * (g$beta / g$se)^2)
    origin <- "Wakefield_ABF_single_signal_sensitivity"
    sd_used <- prior_sd
  }
  if (any(!is.finite(bf))) stop("Invalid GWAS ABF")
  lbf <- matrix(bf, 1, dimnames = list("single_signal", g$key))
  list(lbf = lbf, alpha = matrix(exp(bf - cg_lse(bf)), 1, dimnames = dimnames(lbf)),
       components = "single_signal", variants = g, analysis = origin, effect_prior_sd = sd_used)
}

# The reference must already be restricted to unrelated, ancestry-matched donors,
# or supply a PLINK keep file. Never infer ancestry from having access to GTEx.
cg_ld_reader <- function(bfile, plink, keep_file = NULL, min_samples = 100L, max_snps = 10000L) {
  if (is.null(bfile) || !nzchar(bfile)) stop("GWAS LD reference is not configured")
  if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam"))))) stop("Missing PLINK reference files")
  bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE, data.table = FALSE)
  names(bim) <- c("chromosome", "id", "cm", "position", "a1", "a2")
  bim$key <- cg_key(bim$chromosome, bim$position, bim$a1, bim$a2)
  bad <- duplicated(bim$key) | duplicated(bim$key, fromLast = TRUE) |
    duplicated(bim$id) | duplicated(bim$id, fromLast = TRUE)
  bim <- bim[!bad, ]
  function(g, dir) {
    b <- bim[match(g$key, bim$key), , drop = FALSE]
    present <- !is.na(b$id)
    g <- g[present, , drop = FALSE]; b <- b[present, , drop = FALSE]
    if (nrow(g) < 100) stop("Fewer than 100 GWAS SNPs in the LD reference")
    if (nrow(g) > max_snps) stop("LD matrix would exceed max_ld_snps; increase explicitly if memory permits")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    ids <- file.path(dir, "variants.txt")
    prefix <- file.path(dir, "reference")
    writeLines(b$id, ids)
    args <- c("--bfile", shQuote(bfile), "--extract", shQuote(ids),
              if (!is.null(keep_file)) c("--keep", shQuote(keep_file)),
              "--export", "A", "--out", shQuote(prefix))
    status <- system2(plink, args, stdout = file.path(dir, "plink_stdout.log"), stderr = file.path(dir, "plink_stderr.log"))
    if (status != 0) stop("PLINK LD export failed; inspect ", dir)
    raw <- data.table::fread(paste0(prefix, ".raw"), data.table = FALSE)
    if (nrow(raw) < min_samples) stop("Too few unrelated reference donors")
    # PLINK export A counts the allele in the column suffix, not necessarily ALT.
    cols1 <- paste(b$id, b$a1, sep = "_"); cols2 <- paste(b$id, b$a2, sep = "_")
    use1 <- cols1 %in% names(raw)
    columns <- ifelse(use1, cols1, cols2)
    counted <- ifelse(use1, b$a1, b$a2)
    if (!all(columns %in% names(raw))) stop("Cannot align PLINK dosage columns")
    X <- as.matrix(raw[, columns, drop = FALSE])
    freq <- colMeans(X, na.rm = TRUE) / 2
    missing_rate <- colMeans(is.na(X))
    flip <- counted != g$effect_allele
    if (any(counted != g$effect_allele & counted != g$other_allele)) stop("LD allele mismatch")
    X[, flip] <- 2 - X[, flip, drop = FALSE]
    freq[flip] <- 1 - freq[flip]
    ok <- is.finite(freq) & freq > .01 & freq < .99 & missing_rate <= .05 &
      (is.na(g$eaf) | abs(g$eaf - freq) <= .2)
    g <- g[ok, , drop = FALSE]; X <- X[, ok, drop = FALSE]
    for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
    ok <- apply(X, 2, sd) > 0
    g <- g[ok, , drop = FALSE]; X <- X[, ok, drop = FALSE]
    colnames(X) <- g$key
    list(gwas = g, R = cor(X), n_reference = nrow(X), n_requested = length(present))
  }
}

cg_susie <- function(g, source, ld, dir, L = 10L, max_s = .1, min_ld_fraction = .8) {
  if (!isTRUE(source$primary_eligible)) stop(source$primary_reason)
  ref <- ld(g, dir)
  d <- ref$gwas
  if (nrow(d) < 100 || nrow(d) / nrow(g) < min_ld_fraction) stop("Insufficient GWAS coverage in LD reference")
  if (!g$key[which.max(abs(g$beta / g$se))] %in% d$key) stop("Top regional GWAS Z-score SNP is absent after LD QC")
  if (any(!is.finite(d$n)) || any(d$n <= 2)) stop("Missing usable GWAS sample size")
  n <- median(d$n)
  if (min(d$n) / max(d$n) < .8) stop("Per-SNP GWAS N varies by >20%; harmonize cohort coverage before RSS")
  z <- d$beta / d$se
  mismatch <- susieR::estimate_s_rss(z, ref$R, n = n)
  if (!is.finite(mismatch) || mismatch > max_s) stop("GWAS/LD inconsistency exceeds max_rss_s: ", signif(mismatch, 3))
  fit <- susieR::susie_rss(z = z, R = ref$R, n = n, L = L,
    estimate_residual_variance = FALSE, max_iter = 1000, coverage = .95, min_abs_corr = .5)
  saveRDS(list(fit = fit, variants = d, n = n, n_reference = ref$n_reference,
               rss_s = mismatch, source = source), file.path(dir, "gwas_fit.rds"))
  if (!isTRUE(fit$converged)) stop("GWAS SuSiE did not converge")
  components <- fit$sets$cs_index
  if (!length(components)) stop("no_gwas_cs: no credible GWAS signal; H0-H4 not fabricated from absent CSs")
  bf <- fit$lbf_variable[components, , drop = FALSE]
  alpha <- fit$alpha[components, , drop = FALSE]
  colnames(bf) <- colnames(alpha) <- d$key
  list(lbf = bf, alpha = alpha, components = components, variants = d,
       analysis = "GWAS_SuSiE_RSS", n = n, rss_s = mismatch, n_reference = ref$n_reference)
}

cg_bundle_keys <- function(x) {
  key <- cg_key(x$variants$chromosome, x$variants$position, x$variants$ref, x$variants$alt)
  if (anyDuplicated(key)) stop("Duplicate biological SNPs in QTL bundle")
  colnames(x$lbf) <- colnames(x$alpha) <- key
  x$variants$key <- key
  x
}

cg_pair <- function(bundle, gw, priors, min_mass = .9, min_snps = 100L) {
  eq <- lapply(bundle[c("em", "vanilla")], cg_bundle_keys)
  common <- Reduce(intersect, list(colnames(eq$em$lbf), colnames(eq$vanilla$lbf), colnames(gw$lbf)))
  if (length(common) < min_snps) stop("Insufficient three-way SNP overlap: ", length(common))
  # Coverage is measured BEFORE intersection, both with native QTL alpha and with
  # the uniform SNP prior used for coloc (native SNP priors can differ).
  mass <- function(x) {
    native <- rowSums(x$alpha[, common, drop = FALSE])
    uniform <- apply(x$lbf, 1, function(b) exp(cg_lse(b[match(common, colnames(x$lbf))]) - cg_lse(b)))
    pmin(native, uniform)
  }
  coverage <- lapply(eq, mass); gm <- mass(gw)
  rows <- variants <- list()
  for (method in names(eq)) {
    q <- eq[[method]]
    for (i in seq_len(nrow(q$lbf))) for (j in seq_len(nrow(gw$lbf))) for (k in seq_len(nrow(priors))) {
      p <- priors[k, ]
      if (length(common) * max(p$p1, p$p2) >= 1) stop("Region too large for chosen per-SNP priors")
      good <- coverage[[method]][i] >= min_mass && gm[j] >= min_mass
      row <- data.frame(method, qtl_component = q$components[i], gwas_component = gw$components[j],
        analysis = gw$analysis, p1 = p$p1, p2 = p$p2, p12 = p$p12, nsnps = length(common),
        qtl_retained_mass = coverage[[method]][i], gwas_retained_mass = gm[j],
        status = if (good) "ok" else "insufficient_posterior_coverage",
        H0 = NA_real_, H1 = NA_real_, H2 = NA_real_, H3 = NA_real_, H4 = NA_real_,
        H3_plus_H4 = NA_real_, H4_given_H3_or_H4 = NA_real_, h4_lead = NA_character_, h4_position = NA_real_)
      cs_info <- bundle[[paste0(method, "_cs")]]
      cs_info <- Filter(function(x) identical(as.integer(x$component), as.integer(q$components[i])), cs_info)
      row$qtl_lead_snp <- if (length(cs_info)) cs_info[[1]]$snp else NA_character_
      row$qtl_lead_coding <- if (length(cs_info)) cs_info[[1]]$coding else NA_character_
      row$em_screen_classification <- NA_character_
      row$min_distance_to_vanilla_leads_bp <- NA_real_
      if (method == "em" && length(cs_info)) {
        distances <- vapply(bundle$vanilla_cs, function(v) if (identical(v$chromosome, cs_info[[1]]$chromosome))
          abs(v$position - cs_info[[1]]$position) else NA_real_, numeric(1))
        row$em_screen_classification <- if (!length(distances)) "no_vanilla_cs" else if (anyNA(distances))
          "invalid_coordinates" else if (all(distances > bundle$metadata$shift_bp))
          "shifted_from_all_vanilla_leads" else "not_shifted"
        if (length(distances) && !anyNA(distances)) row$min_distance_to_vanilla_leads_bp <- min(distances)
      }
      if (good) {
        b1 <- setNames(q$lbf[i, common], common); b2 <- setNames(gw$lbf[j, common], common)
        result <- coloc::coloc.bf_bf(b1, b2, p1 = p$p1, p2 = p$p2, p12 = p$p12,
                                    trim_by_posterior = FALSE)
        pp <- unlist(as.data.frame(result$summary)[1, paste0("PP.H", 0:4, ".abf")], use.names = FALSE)
        if (any(!is.finite(pp)) || abs(sum(pp) - 1) > 1e-8) stop("Invalid coloc posterior")
        row[paste0("H", 0:4)] <- as.list(pp)
        row$H3_plus_H4 <- pp[4] + pp[5]
        row$H4_given_H3_or_H4 <- if (pp[4] + pp[5] >= .1) pp[5] / (pp[4] + pp[5]) else NA_real_
        joint <- b1 + b2
        snp_pp <- exp(joint - cg_lse(joint))
        lead <- names(which.max(snp_pp))
        row$h4_lead <- lead
        row$h4_position <- gw$variants$position[match(lead, gw$variants$key)]
        variants[[length(variants) + 1L]] <- data.frame(method, qtl_component = q$components[i],
          gwas_component = gw$components[j], p1 = p$p1, p2 = p$p2, p12 = p$p12,
          snp = common, SNP_PP_H4 = unname(snp_pp))
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  list(pairs = data.table::rbindlist(rows, fill = TRUE),
       variants = data.table::rbindlist(variants, fill = TRUE))
}

cg_changes <- function(pairs, h4_threshold = .8, min_delta = .1) {
  if (!nrow(pairs)) return(data.frame())
  groups <- unique(pairs[, c("gwas_component", "analysis", "p1", "p2", "p12"), with = FALSE])
  rows <- list()
  for (k in seq_len(nrow(groups))) {
    group <- groups[k, ]
    d <- pairs[gwas_component == group$gwas_component & analysis == group$analysis &
                 p1 == group$p1 & p2 == group$p2 & p12 == group$p12]
    complete <- all(c("em", "vanilla") %in% d$method) && all(d$status == "ok")
    row <- as.data.frame(group)
    row$summary_rule <- "max_H4_pair_per_method_for_same_GWAS_signal"
    row$status <- if (complete) "ok" else "unassessed_missing_or_filtered_component"
    row$classification <- "unassessed"
    if (complete) {
      e <- d[method == "em"]; e <- e[which.max(H4)]
      v <- d[method == "vanilla"]; v <- v[which.max(H4)]
      for (h in paste0("H", 0:4)) {
        row[[paste0("em_", h)]] <- e[[h]]; row[[paste0("vanilla_", h)]] <- v[[h]]
        row[[paste0("delta_", h)]] <- e[[h]] - v[[h]]
      }
      row$em_component <- e$qtl_component; row$vanilla_component <- v$qtl_component
      row$em_screen_classification <- e$em_screen_classification
      row$em_lead_coding <- e$qtl_lead_coding
      row$em_h4_lead <- e$h4_lead; row$vanilla_h4_lead <- v$h4_lead
      row$h4_lead_distance_bp <- abs(e$h4_position - v$h4_position)
      row$vanilla_top_hypothesis <- paste0("H", which.max(unlist(v[, paste0("H", 0:4), with = FALSE])) - 1L)
      row$em_top_hypothesis <- paste0("H", which.max(unlist(e[, paste0("H", 0:4), with = FALSE])) - 1L)
      row$classification <- if (e$H4 >= h4_threshold && v$H4 < h4_threshold && row$delta_H4 >= min_delta) "gain" else
        if (v$H4 >= h4_threshold && e$H4 < h4_threshold && row$delta_H4 <= -min_delta) "loss" else
        if (e$H4 >= h4_threshold && v$H4 >= h4_threshold) "both_supported" else "neither_crosses_threshold"
    }
    rows[[length(rows) + 1L]] <- row
  }
  as.data.frame(data.table::rbindlist(rows, fill = TRUE))
}

run_em_gwas_coloc <- function(job_dir, modes = c("abf", "susie"), traits = NULL,
    ld_bfile = NULL, ld_keep_file = NULL, ld_ancestry = NULL, ld_build = NULL,
    plink = "plink2", min_reference_samples = 100L, max_ld_snps = 10000L,
    L = 10L, max_rss_s = .1, min_ld_fraction = .8, min_mass = .9, min_snps = 100L,
    min_info = .8, priors = data.frame(p1 = 1e-4, p2 = 1e-4, p12 = c(1e-6, 5e-6, 1e-5)),
    h4_threshold = .8, min_delta = .1, effect_prior_sd = list()) {
  cg_need(c("data.table", "jsonlite", "digest", "R.utils", "coloc"))
  if (any(!modes %in% c("abf", "susie"))) stop("modes must be abf and/or susie")
  if (any(!is.finite(as.matrix(priors))) || any(priors <= 0) || any(priors$p12 > pmin(priors$p1, priors$p2))) stop("Invalid coloc priors")
  job <- readRDS(file.path(job_dir, "job_metadata.rds"))
  cases <- read.csv(file.path(job_dir, "cases.csv"), stringsAsFactors = FALSE)
  regions <- read.csv(file.path(job_dir, "regions.csv"), stringsAsFactors = FALSE)
  sources <- jsonlite::fromJSON(file.path(job_dir, "sources.json"), simplifyVector = FALSE)$sources
  if (!is.null(traits)) {
    if (any(!traits %in% vapply(sources, `[[`, "", "trait_id"))) stop("Unknown trait ID")
    sources <- Filter(function(s) s$trait_id %in% traits, sources)
  }
  out <- tempfile(paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = job_dir)
  dir.create(out)
  settings <- as.list(environment()); settings <- settings[setdiff(names(formals(sys.function())), c("job_dir"))]
  saveRDS(list(settings = settings, job = job, session = sessionInfo()), file.path(out, "analysis_metadata.rds"))
  ld <- NULL; ld_problem <- "LD not configured; set ancestry, GRCh38 PLINK prefix and donor subset"
  if ("susie" %in% modes && !is.null(ld_bfile)) {
    cg_need("susieR")
    if (!identical(ld_build, "GRCh38") || is.null(ld_ancestry)) stop("Declare LD build GRCh38 and ancestry explicitly")
    ld <- cg_ld_reader(ld_bfile, plink, ld_keep_file, min_reference_samples, max_ld_snps)
  }
  all_pairs <- changes <- audit <- list()
  record <- function(trait, case, mode, status, detail) {
    audit[[length(audit) + 1L]] <<- data.frame(trait_id = trait, case_id = case, mode, status, detail)
  }
  for (s in sources) {
    trait <- s$trait_id
    message("Coloc: ", trait)
    file <- file.path(job_dir, "gwas", paste0(trait, ".tsv.gz"))
    input <- tryCatch({
      if (!file.exists(file) || !file.exists(file.path(job_dir, "gwas", paste0(trait, ".json"))))
        stop("Download/extraction is not complete for ", trait, "; inspect gwas/download_status.json")
      receipt <- jsonlite::fromJSON(file.path(job_dir, "gwas", paste0(trait, ".json")), simplifyVector = FALSE)
      if (receipt$output_sha256 != digest::digest(file = file, algo = "sha256")) stop("GWAS regional checksum mismatch")
      if (!identical(s, receipt$source)) stop("Frozen source definition disagrees with downloaded data")
      data.table::fread(file, data.table = FALSE)
    }, error = identity)
    for (i in seq_len(nrow(cases))) {
      case <- cases$case_id[i]
      prefix <- list(trait_id = trait, case_id = case, gene = cases$gene[i], tissue = cases$tissue[i], em_iteration = job$iteration)
      if (inherits(input, "error")) {
        for (mode in modes) record(trait, case, mode, "missing_or_invalid_gwas", conditionMessage(input))
        next
      }
      region <- regions[match(case, regions$case_id), ]
      raw <- input[input$chromosome == region$chromosome & input$position >= region$start & input$position <= region$end, , drop = FALSE]
      cleaned <- cg_clean_gwas(raw, min_info)
      g <- cleaned$data
      bundle_path <- file.path(job_dir, cases$bundle[i])
      if (unname(tools::md5sum(bundle_path)) != cases$md5[i]) stop("Frozen BF bundle checksum mismatch")
      bundle <- readRDS(bundle_path)
      for (mode in modes) {
        if (mode == "susie" && (is.null(ld) || !identical(s$ancestry, ld_ancestry))) {
          record(trait, case, mode, "pending_ld", if (is.null(ld)) ld_problem else "GWAS and reference ancestry disagree")
          next
        }
        dir <- file.path(out, trait, case, mode)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        result <- tryCatch({
          if (nrow(g) < min_snps) stop("Insufficient regional GWAS SNPs after allele/INFO QC")
          gw <- if (mode == "abf") cg_abf(g, s, effect_prior_sd[[trait]]) else
            cg_susie(g, s, ld, dir, L, max_rss_s, min_ld_fraction)
          saveRDS(gw, file.path(dir, "gwas_bfs.rds"))
          paired <- cg_pair(bundle, gw, priors, min_mass, min_snps)
          pairs <- cbind(as.data.frame(prefix), as.data.frame(paired$pairs))
          delta <- cg_changes(paired$pairs, h4_threshold, min_delta)
          if (nrow(delta)) delta <- cbind(as.data.frame(prefix), delta)
          data.table::fwrite(paired$variants, file.path(dir, "snp_h4_posteriors.tsv.gz"), sep = "\t")
          saveRDS(list(source = s, qc = cleaned$counts, gwas_n = gw$n, rss_s = gw$rss_s,
            effect_prior_sd = gw$effect_prior_sd), file.path(dir, "qc.rds"))
          list(pairs = pairs, changes = delta)
        }, error = identity)
        if (inherits(result, "error")) record(trait, case, mode, "unassessed", conditionMessage(result)) else {
          all_pairs[[length(all_pairs) + 1L]] <- result$pairs
          if (nrow(result$changes)) changes[[length(changes) + 1L]] <- result$changes
          record(trait, case, mode, "completed", paste(nrow(result$pairs), "signal/prior comparisons"))
        }
      }
    }
    cg_write(all_pairs, file.path(out, "all_signal_pairs.csv"))
    cg_write(changes, file.path(out, "paired_changes.csv"))
    cg_write(audit, file.path(out, "audit.csv"))
  }
  ch <- cg_write(changes, file.path(out, "paired_changes.csv"))
  cg_write(all_pairs, file.path(out, "all_signal_pairs.csv"))
  cg_write(audit, file.path(out, "audit.csv"))
  if (nrow(ch)) {
    counts <- ch[, .N, by = .(trait_id, analysis, p1, p2, p12, classification)]
    data.table::fwrite(counts, file.path(out, "change_counts.csv"))
  }
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(out, "completed.txt"))
  message("Results and all unassessed cases: ", out)
  invisible(out)
}
