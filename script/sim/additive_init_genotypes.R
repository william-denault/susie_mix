# Use the original .raw pool when present; otherwise extract a real GTEx gene
# region with the scanning workhorse's paths, annotation parser and PLINK flags.
init_file_signatures <- function(files) {
  info <- file.info(files)
  data.frame(file = normalizePath(files, winslash = "/", mustWork = TRUE),
             size = info$size, mtime = as.numeric(info$mtime), row.names = NULL)
}

prepare_additive_genotypes <- function(
    genotype_dir = Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(init_project_dir, "temp_plink")),
    datadir = Sys.getenv("SUSIE_MIX_GTEX_DIR", "/project2/mstephens/gtex"),
    project_dir = init_project_dir,
    target_gene = Sys.getenv("SUSIE_MIX_SIM_GENE", ""),
    cis_window = 5e5, plink_threads = 1L,
    plink_memory = as.numeric(Sys.getenv("SUSIE_MIX_PLINK_MEMORY_MB", "8000"))) {
  if (!nzchar(genotype_dir)) genotype_dir <- file.path(init_project_dir, "temp_plink")
  files <- list.files(genotype_dir, pattern = "\\.raw$", full.names = TRUE)
  if (length(files)) {
    if (any(file.access(files, 4L) != 0L)) stop("Some simulation genotype files are not readable: ", genotype_dir)
    directory <- normalizePath(genotype_dir, winslash = "/", mustWork = TRUE)
    files <- normalizePath(files, winslash = "/", mustWork = TRUE)
    message("Using the original simulation genotype pool: ", length(files), " regions in ", directory)
    return(list(mode = "raw", directory = directory, files = files,
         settings = list(mode = "original_simulation_files", files = init_file_signatures(files))))
  }

  stopifnot(length(target_gene) == 1L, !is.na(target_gene),
    cis_window > 0, cis_window == floor(cis_window),
    plink_threads >= 1L, plink_threads == floor(plink_threads),
    length(plink_memory) == 1L, is.finite(plink_memory),
    plink_memory >= 640L, plink_memory <= .Machine$integer.max, plink_memory == floor(plink_memory))
  bfile <- file.path(datadir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv")
  gtf <- file.path(datadir, "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
  plink <- Sys.getenv("SUSIE_MIX_PLINK", file.path(datadir, "plink2"))
  gene_file <- file.path(project_dir, "data/genes_protein_coding.txt")
  annotation_script <- file.path(project_dir, "script/scan_tissue_attempt/get_gene_annotations.R")
  required <- c(paste0(bfile, c(".bed", ".bim", ".fam")), gtf, plink, annotation_script,
                if (!nzchar(target_gene)) gene_file)
  missing <- required[!file.exists(required) | file.access(required, 4L) != 0L]
  if (length(missing)) stop("No cached .raw files; GTEx extraction needs these missing/unreadable files:\n",
    paste(missing, collapse = "\n"), "\nThe default data directory matches the scanning workhorse: ", datadir,
    call. = FALSE)
  if (.Platform$OS.type != "windows" && file.access(plink, 1L) != 0L)
    stop("PLINK is not executable: ", plink)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install required package: data.table")
  annotation_env <- new.env(parent = environment())
  annotation_env$fread <- data.table::fread
  source(annotation_script, local = annotation_env)
  annotations <- annotation_env$get_gene_annotations(gtf)
  genes <- if (nzchar(target_gene)) target_gene else unique(trimws(readLines(gene_file, warn = FALSE)))
  genes <- genes[nzchar(genes)]
  # As in run_susie_gene(), use the first matching annotation for a gene.
  annotations <- annotations[match(genes, annotations$gene_name), , drop = FALSE]
  chr <- sub("^chr", "", annotations$chromosome)
  keep <- !is.na(annotations$gene_name) & chr %in% as.character(1:22) &
    annotations$strand %in% c("+", "-") & is.finite(annotations$start) & is.finite(annotations$end)
  annotations <- annotations[keep, , drop = FALSE]
  chr <- chr[keep]
  tss <- ifelse(annotations$strand == "+", annotations$start, annotations$end)
  loci <- data.frame(gene = annotations$gene_name, chromosome = chr, tss = tss,
    from_bp = pmax(0, tss - cis_window), to_bp = tss + cis_window, row.names = NULL)
  if (!nrow(loci)) stop("No matching autosomal protein-coding gene regions in the GTEx annotations.")
  message("No cached .raw files: using real GTEx cis-genotypes from ", datadir, ". ",
    if (nzchar(target_gene)) paste0("Gene: ", target_gene) else
      paste0("Each seed selects one of ", nrow(loci), " autosomal protein-coding genes."),
    " PLINK workspace: ", plink_memory, " MiB.")
  list(mode = "gtex", bfile = bfile, plink = plink, loci = loci,
    plink_threads = as.integer(plink_threads), plink_memory = as.integer(plink_memory),
    settings = list(mode = "gtex_fallback", files = init_file_signatures(required), loci = loci,
      cis_window = cis_window, min_maf_plink = 0, target_gene = target_gene,
      plink_threads = as.integer(plink_threads), plink_memory = as.integer(plink_memory)))
}

.init_plink_call <- function(executable, args, logfile) {
  system2(executable, args = vapply(args, shQuote, ""), stdout = logfile, stderr = logfile)
}

extract_additive_region <- function(source, index, temp_dir = tempdir(), run_plink = .init_plink_call) {
  stopifnot(identical(source$mode, "gtex"), length(index) == 1L,
    index >= 1L, index <= nrow(source$loci), index == floor(index))
  locus <- source$loci[index, , drop = FALSE]
  # A private folder lets the original sim_mix() read just this gene's export.
  directory <- tempfile("additive_init_gene_", tmpdir = temp_dir)
  if (!dir.create(directory, recursive = TRUE)) stop("Cannot create genotype extraction directory: ", directory)
  prefix <- file.path(directory, "genotypes")
  cleanup <- function() {
    # Only this extraction's files; R removes the empty temporary folder on exit.
    paths <- Sys.glob(paste0(prefix, ".*"))
    if (length(paths)) unlink(paths, recursive = FALSE)
    invisible(NULL)
  }
  succeeded <- FALSE
  on.exit(if (!succeeded) cleanup(), add = TRUE)
  args <- c("--bfile", source$bfile, "--chr", locus$chromosome,
    "--from-bp", format(locus$from_bp, scientific = FALSE, trim = TRUE),
    "--to-bp", format(locus$to_bp, scientific = FALSE, trim = TRUE),
    "--snps-only", "--max-alleles", "2", "--rm-dup", "exclude-all",
    "--threads", as.character(source$plink_threads), "--memory", as.character(source$plink_memory),
    "--maf", "0", "--recode", "A", "--out", prefix)
  logfile <- paste0(prefix, ".console.log")
  message("Extracting GTEx gene ", locus$gene, " (chr", locus$chromosome, ":",
          locus$from_bp, "-", locus$to_bp, ")")
  status <- run_plink(source$plink, args, logfile)
  raw_file <- paste0(prefix, ".raw")
  if (status != 0L || !file.exists(raw_file)) {
    log_lines <- if (file.exists(logfile)) readLines(logfile, warn = FALSE) else character()
    detail <- paste(tail(log_lines, 12L), collapse = "\n")
    if (any(grepl("out of memory|cannot allocate|failed to allocate|bad_alloc", log_lines, ignore.case = TRUE))) {
      message <- paste0("PLINK ran out of memory for ", locus$gene,
        " with --memory ", source$plink_memory, " MiB (exit ", status, ").\n",
        "Stopping the experiment. Increase SUSIE_MIX_PLINK_MEMORY_MB within the RStudio/job memory allocation, then rerun.\n",
        detail)
      stop(structure(list(message = message, call = NULL),
        class = c("init_plink_memory_error", "error", "condition")))
    }
    stop("PLINK extraction failed for ", locus$gene, " (exit ", status, ").\n", detail)
  }
  succeeded <- TRUE
  list(file = raw_file, directory = directory, gene = locus$gene, region = locus, cleanup = cleanup)
}
