# Read-only against source data and EM outputs. Reports are confined to bundle/checks.
stable_id <- function(x) sub("\\.[0-9]+$", "", x)
read_rename_map <- function(bundle) {
  x <- read.csv(file.path(bundle, "renamed_gene_map.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(x) == 680L, length(unique(x$expression_symbol)) == 679L,
            !anyDuplicated(x$gene), !anyDuplicated(x$stable_id),
            all(grepl("^ENSG[0-9]+$", x$stable_id)),
            all(x$chromosome %in% as.character(1:22)))
  x
}

check_metadata <- function(map, expression, annotation) {
  expression_ids <- stable_id(expression$Name)
  annotation_ids <- stable_id(annotation$ensembl)
  checks <- map
  checks$expression_matches <- vapply(map$stable_id, function(id) sum(expression_ids == id), integer(1))
  checks$annotation_matches <- vapply(map$stable_id, function(id) sum(annotation_ids == id), integer(1))
  checks$actual_expression_id <- checks$actual_expression_symbol <- ""
  checks$actual_annotation_symbol <- checks$actual_chromosome <- ""
  checks$actual_start <- checks$actual_end <- NA_real_
  checks$issue <- ""
  for (i in seq_len(nrow(map))) {
    issue <- character()
    e <- which(expression_ids == map$stable_id[i]); a <- which(annotation_ids == map$stable_id[i])
    if (length(e) != 1L) issue <- c(issue, "Expression ID missing or ambiguous") else {
      checks$actual_expression_id[i] <- expression$Name[e]
      checks$actual_expression_symbol[i] <- expression$Description[e]
      if (!identical(expression$Name[e], map$expression_id[i]) ||
          !identical(expression$Description[e], map$expression_symbol[i]))
        issue <- c(issue, "Expression ID version or symbol differs from audited GTEx v8 row")
    }
    if (length(a) != 1L) issue <- c(issue, "Protein-coding annotation ID missing or ambiguous") else {
      chr <- sub("^chr", "", as.character(annotation$chromosome[a]))
      checks$actual_annotation_symbol[i] <- annotation$gene_name[a]
      checks$actual_chromosome[i] <- chr
      checks$actual_start[i] <- annotation$start[a]
      checks$actual_end[i] <- annotation$end[a]
      if (chr != map$chromosome[i] || annotation$start[a] != map$start[i] ||
          annotation$end[a] != map$end[i] || annotation$strand[a] != map$strand[i] ||
          annotation$gene_name[a] != map$annotation_symbol[i])
        issue <- c(issue, "RCC annotation differs from the audited Ensembl 103 locus")
    }
    checks$issue[i] <- paste(issue, collapse = "; ")
  }
  checks$status <- ifelse(nzchar(checks$issue), "FAIL", "PASS")
  checks
}

run_preflight <- function(bundle, project, datadir = "/project2/mstephens/gtex") {
  bundle <- normalizePath(bundle, winslash = "/", mustWork = TRUE)
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  map <- read_rename_map(bundle)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  expression_file <- file.path(datadir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
  annotation_file <- file.path(datadir, "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
  cat("Reading expression row IDs/names only. No expression values or model fitting.\n")
  expression <- as.data.frame(data.table::fread(expression_file, skip = 2, select = c("Name", "Description")))
  helper <- new.env(parent = globalenv()); helper$fread <- data.table::fread
  sys.source(file.path(bundle, "script/scan_tissue_attempt/get_gene_annotations.R"), helper)
  annotation <- helper$get_gene_annotations(annotation_file)
  checks <- check_metadata(map, expression, annotation)
  manifest_files <- c(mix = file.path(project, "results_em/iteration_018/manifest.csv"),
                      slider = file.path(project, "results_slide_prior_em/iteration_000/manifest.csv"))
  if (any(!file.exists(manifest_files))) stop("Original mix/slider manifest missing.")
  manifests <- lapply(manifest_files, function(p) read.csv(p, stringsAsFactors = FALSE)$gene)
  if (!setequal(manifests$mix, manifests$slider)) stop("Original mix and slider gene cohorts differ.")
  if (length(intersect(map$expression_symbol, manifests$mix))) stop("A recovery symbol is already in the old cohort.")
  if (anyDuplicated(manifests$mix)) stop("Duplicate original manifest gene.")
  dirs <- list.files(file.path(project, "results_slide_prior_em"), "^iteration_[0-9]+$", full.names = TRUE)
  completed <- dirs[file.exists(file.path(dirs, "recovery/parallel/COMPLETE.rds"))]
  completed <- completed[order(as.integer(sub("iteration_", "", basename(completed))))]
  checkpoints <- data.frame(model = c("mix", "slider"),
     path = c(file.path(project, "results_em/iteration_018"), if (length(completed)) tail(completed, 1) else ""))
  checks_dir <- file.path(bundle, "checks"); dir.create(checks_dir, showWarnings = FALSE)
  write.csv(checks, file.path(checks_dir, "mapping_checks.csv"), row.names = FALSE)
  write.csv(checkpoints, file.path(checks_dir, "pilot_checkpoints.csv"), row.names = FALSE)
  stat_files <- c(expression_file, annotation_file, manifest_files)
  input_info <- file.info(stat_files)[, c("size", "mtime"), drop = FALSE]
  input_info$path <- rownames(input_info); rownames(input_info) <- NULL
  packages <- vapply(c("susieR", "susieRSlidePrior", "data.table", "matrixStats"), function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "MISSING", character(1))
  state <- list(project = project, datadir = datadir, map_md5 = unname(tools::md5sum(file.path(bundle, "renamed_gene_map.csv"))),
                annotation_md5 = unname(tools::md5sum(annotation_file)), input_info = input_info,
                packages = packages, checked_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                passed = all(checks$status == "PASS"))
  saveRDS(state, file.path(checks_dir, "preflight.rds"))
  lines <- c(sprintf("Candidate symbols: %d", length(unique(map$expression_symbol))),
      sprintf("Candidate expression IDs: %d", nrow(map)),
      sprintf("Exact expression and annotation matches: %d", sum(checks$status == "PASS")),
      sprintf("Failed/ambiguous mappings: %d", sum(checks$status == "FAIL")),
      "C2orf61 has two distinct GTEx IDs, retained with separate filenames.",
      "CTAGE5 is selected by its exact GTEx ID; existing MIA2 results are never overwritten.",
      "Input presence does not guarantee eligible tissue fits or convergence.",
      paste(names(packages), packages, sep = ": "))
  writeLines(lines, file.path(checks_dir, "summary.txt")); cat(paste(lines, collapse = "\n"), "\n")
  if (!state$passed) stop("Preflight failed. Inspect mapping_checks.csv; no EM files changed.")
  invisible(checks)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) < 2 || length(args) > 3) stop("Usage: Rscript preflight.R BUNDLE PROJECT [GTEX_DATA_DIR]")
  run_preflight(args[1], args[2], if (length(args)==3) args[3] else "/project2/mstephens/gtex")
}
