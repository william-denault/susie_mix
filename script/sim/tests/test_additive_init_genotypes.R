source("script/sim/sim_additive_slide_init.R")
test_original_genotype_pool <- function() {
  directory <- tempfile("original_genotype_pool_")
  dir.create(directory)
  failure <- tryCatch(prepare_additive_genotypes(directory, datadir = file.path(directory, "missing_GTEx")),
                      error = conditionMessage)
  stopifnot(grepl("GTEx extraction needs", failure, fixed = TRUE))
  files <- file.path(directory, c("b.raw", "a.raw"))
  for (file in files) writeLines("fixture", file)
  source <- prepare_additive_genotypes(directory)
  stopifnot(source$mode == "raw",
    identical(source$files, normalizePath(sort(files), winslash = "/", mustWork = TRUE)),
    source$settings$mode == "original_simulation_files")
  # Reject a changing locus pool before attempting to fit any data.
  writeLines("new region", file.path(directory, "c.raw"))
  failure <- tryCatch(simulate_additive_init_data(1L, genotype_source = source), error = conditionMessage)
  stopifnot(grepl("file list changed", failure, fixed = TRUE))
  cat("PASS: original genotype pool takes precedence, missing GTEx diagnostics, and changed-pool detection.\n")
}
test_original_genotype_pool()

test_gt_extraction <- function() {
  directory <- tempfile("additive_gtex_validation_")
  datadir <- file.path(directory, "GTEx data")
  project <- file.path(directory, "project")
  for (path in c(datadir, file.path(project, "data"),
                 file.path(project, "script/scan_tissue_attempt")))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  bfile <- file.path(datadir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv")
  for (path in c(paste0(bfile, c(".bed", ".bim", ".fam")), file.path(datadir, "plink2")))
    writeLines("fixture", path)
  Sys.chmod(file.path(datadir, "plink2"), "0755")
  writeLines(c("A", "B", "C", "Missing"), file.path(project, "data/genes_protein_coding.txt"))
  file.copy(file.path(init_project_dir, "script/scan_tissue_attempt/get_gene_annotations.R"),
            file.path(project, "script/scan_tissue_attempt/get_gene_annotations.R"))
  gtf <- file.path(datadir, "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
  connection <- gzfile(gtf, "wt")
  writeLines(c("# fixture",
    'chr1\tensembl\tgene\t1000000\t1010000\t.\t+\t.\tgene_id "ENSA"; gene_type "protein_coding"; gene_name "A";',
    'chr2\tensembl\tgene\t2000000\t2010000\t.\t-\t.\tgene_id "ENSB"; gene_type "protein_coding"; gene_name "B";',
    'chrX\tensembl\tgene\t3000000\t3010000\t.\t+\t.\tgene_id "ENSC"; gene_type "protein_coding"; gene_name "C";'), connection)
  close(connection)
  previous_plink <- Sys.getenv("SUSIE_MIX_PLINK", unset = NA_character_)
  Sys.unsetenv("SUSIE_MIX_PLINK")
  on.exit(if (is.na(previous_plink)) Sys.unsetenv("SUSIE_MIX_PLINK") else
    Sys.setenv(SUSIE_MIX_PLINK = previous_plink), add = TRUE)
  previous_memory <- Sys.getenv("SUSIE_MIX_PLINK_MEMORY_MB", unset = NA_character_)
  Sys.unsetenv("SUSIE_MIX_PLINK_MEMORY_MB")
  on.exit(if (is.na(previous_memory)) Sys.unsetenv("SUSIE_MIX_PLINK_MEMORY_MB") else
    Sys.setenv(SUSIE_MIX_PLINK_MEMORY_MB = previous_memory), add = TRUE)
  prepare <- function(target_gene = "") prepare_additive_genotypes(
    genotype_dir = file.path(directory, "empty"), datadir = datadir,
    project_dir = project, target_gene = target_gene)
  source <- prepare()
  stopifnot(identical(source$mode, "gtex"), identical(source$loci$gene, c("A", "B")),
    identical(source$loci$chromosome, c("1", "2")), source$plink_threads == 1L,
    source$plink_memory == 8000L,
    identical(as.numeric(source$loci$tss), c(1000000, 2010000)),
    identical(as.numeric(source$loci$from_bp), c(500000, 1510000)),
    identical(prepare("B")$loci$gene, "B"))
  failure <- tryCatch(prepare("Missing"), error = conditionMessage)
  stopifnot(grepl("No matching autosomal", failure, fixed = TRUE))
  Sys.setenv(SUSIE_MIX_PLINK_MEMORY_MB = "16000")
  stopifnot(prepare()$plink_memory == 16000L)
  Sys.unsetenv("SUSIE_MIX_PLINK_MEMORY_MB")

  genotypes <- file.path(directory, "genotypes")
  dir.create(genotypes)
  set.seed(9014)
  X <- sapply(seq(.08, .45, length.out = 24), function(maf) rbinom(650, 2, maf))
  colnames(X) <- paste0("snp", seq_len(ncol(X)))
  raw <- data.frame(FID = 1:650, IID = paste0("donor", 1:650), PAT = 0, MAT = 0,
                    SEX = 0, PHENOTYPE = -9, X)
  raw_file <- file.path(genotypes, "fixture.raw")
  write.table(raw, raw_file, row.names = FALSE, quote = FALSE, sep = "\t")
  called <- NULL
  produced <- character()
  fake_plink <- function(executable, args, logfile) {
    called <<- args
    prefix <- args[match("--out", args) + 1L]
    produced <<- c(produced, paste0(prefix, ".raw"))
    file.copy(raw_file, paste0(prefix, ".raw"))
    writeLines("mock extraction succeeded", logfile)
    0L
  }
  unrelated <- file.path(directory, "unrelated.raw")
  writeLines("keep", unrelated)
  region <- extract_additive_region(source, 2L, temp_dir = directory, run_plink = fake_plink)
  value <- function(flag) called[match(flag, called) + 1L]
  stopifnot(file.exists(region$file), region$gene == "B", value("--bfile") == bfile,
    value("--chr") == "2", value("--from-bp") == "1510000", value("--to-bp") == "2510000",
    value("--recode") == "A", value("--maf") == "0", value("--threads") == "1",
    value("--memory") == "8000")
  region$cleanup()
  stopifnot(!file.exists(region$file), !length(list.files(region$directory)), file.exists(unrelated))
  failed_prefix <- NULL
  failure <- tryCatch(extract_additive_region(source, 1L, temp_dir = directory,
    run_plink = function(executable, args, logfile) {
      failed_prefix <<- args[match("--out", args) + 1L]
      writeLines("partial output", paste0(failed_prefix, ".raw"))
      writeLines("test PLINK failure detail", logfile)
      7L
    }), error = conditionMessage)
  stopifnot(grepl("exit 7", failure, fixed = TRUE),
    grepl("test PLINK failure detail", failure, fixed = TRUE),
    !length(Sys.glob(paste0(failed_prefix, ".*"))), file.exists(unrelated))

  original_call <- .init_plink_call
  assign(".init_plink_call", fake_plink, envir = .GlobalEnv)
  on.exit(assign(".init_plink_call", original_call, envir = .GlobalEnv), add = TRUE)
  direct <- simulate_additive_init_data(1000001L, pve = .025, genotype_source = source)
  cached <- simulate_additive_init_data(1000001L, genotype_dir = genotypes, pve = .025)
  higher_pve <- simulate_additive_init_data(1000001L, pve = .05, genotype_source = source)
  repeated <- simulate_additive_init_data(1000001L, pve = .025, genotype_source = source)
  stopifnot(identical(direct$X, cached$X), identical(direct$y, cached$y),
    identical(direct$true_pos, cached$true_pos), direct$genotype_mode == "gtex_fallback",
    identical(direct$X, higher_pve$X), identical(direct$true_pos, higher_pve$true_pos),
    identical(direct$gene, higher_pve$gene), identical(direct$gene, repeated$gene),
    identical(direct$y, repeated$y), !any(file.exists(produced)), file.exists(cached$raw_file))
  noise <- function(x, pve) drop(x$y - scale(x$X[, x$true_pos]) %*% x$beta_standardized) / sqrt(1 - pve)
  stopifnot(isTRUE(all.equal(noise(direct, .025), noise(higher_pve, .05))))
  # QC/data-generation errors also release the extraction files.
  failure <- tryCatch(simulate_additive_init_data(1000001L, n = 800L, genotype_source = source),
                      error = conditionMessage)
  stopifnot(grepl("Not enough donors", failure, fixed = TRUE), !any(file.exists(produced)))
  cat("PASS: GTEx fallback, optional gene, workhorse cis windows/flags, cleanup, unchanged sim_mix generation, paired PVEs and reproducibility.\n")

  # Reproduce the reported 2000-MiB failure. It must stop the whole experiment
  # after one extraction, checkpoint the failed seed, and release its files.
  source$plink_memory <- source$settings$plink_memory <- 2000L
  original_prepare <- prepare_additive_genotypes
  assign("prepare_additive_genotypes", function(...) source, envir = .GlobalEnv)
  on.exit(assign("prepare_additive_genotypes", original_prepare, envir = .GlobalEnv), add = TRUE)
  extraction_attempts <- 0L
  assign(".init_plink_call", function(executable, args, logfile) {
    extraction_attempts <<- extraction_attempts + 1L
    prefix <- args[match("--out", args) + 1L]
    produced <<- c(produced, paste0(prefix, ".raw"))
    writeLines("partial", paste0(prefix, ".raw"))
    writeLines(c("reserving 2000 MiB for main workspace.", "Error: Out of memory."), logfile)
    2L
  }, envir = .GlobalEnv)
  experiment_dir <- file.path(directory, "oom_resume")
  failure <- tryCatch(run_additive_initialization_experiment(chunks = 1:2,
    reps_per_chunk = 2L, output_dir = experiment_dir), error = identity)
  checkpoint <- file.path(experiment_dir, "pve_0.025/chunks/additive_init_chunk001.rds")
  saved <- readRDS(checkpoint)
  stopifnot(inherits(failure, "init_plink_memory_error"), extraction_attempts == 1L,
    !any(file.exists(produced)), !dir.exists(file.path(experiment_dir, "pve_0.05")),
    !file.exists(file.path(dirname(checkpoint), "additive_init_chunk002.rds")),
    grepl("2000", saved$results[[1L]]$error, fixed = TRUE), is.null(saved$results[[2L]]))
  # Include an old-format failed record, as produced before fail-fast handling.
  saved$results[[2L]] <- list(seed = 1000002L, error = "PLINK extraction failed: Out of memory.")
  saveRDS(saved, checkpoint)
  assign(".init_plink_call", fake_plink, envir = .GlobalEnv)
  source$plink_memory <- source$settings$plink_memory <- 8000L
  output_dir <- file.path(experiment_dir, "pve_0.025")
  run_additive_initialization(chunk = 1L, reps_per_chunk = 2L, output_dir = output_dir,
    pve = .025, genotype_source = source)
  recovered <- readRDS(checkpoint)
  stopifnot(recovered$settings$genotype_inputs$plink_memory == 8000L,
    all(vapply(recovered$results, function(x) is.null(x$error), logical(1))))
  # A further memory-only change skips successes; mixed-memory chunks summarize.
  digest <- tools::md5sum(checkpoint)
  source$plink_memory <- source$settings$plink_memory <- 16000L
  run_additive_initialization(chunk = 1L, reps_per_chunk = 2L, output_dir = output_dir,
    pve = .025, genotype_source = source)
  stopifnot(identical(digest, tools::md5sum(checkpoint)))
  run_additive_initialization(chunk = 2L, reps_per_chunk = 2L, output_dir = output_dir,
    pve = .025, genotype_source = source)
  summary <- summarize_additive_initialization(output_dir)
  stopifnot(nrow(summary$metrics) == 12L, nrow(summary$errors) == 0L)
  cat("PASS: memory override, immediate OOM stop, saved failure, retry of old 2000-MiB checkpoints, preserved successes and mixed-memory summaries.\n")
}
test_gt_extraction()
