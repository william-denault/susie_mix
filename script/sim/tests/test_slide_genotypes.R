# Run from the project root. Uses real five-method fits and a mocked PLINK
# export, because the protected GTEx data are only available on RCC.
source("script/sim/run_job.R")

test_slide_genotypes <- function() {
  directory <- tempfile("slide_gtex_")
  project <- file.path(directory, "project")
  datadir <- file.path(directory, "GTEx data")
  raw_dir <- file.path(directory, "cached genotypes")
  for (path in c(datadir, raw_dir, file.path(project, "data"),
                 file.path(project, "script/scan_tissue_attempt")))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy("script/scan_tissue_attempt/get_gene_annotations.R",
    file.path(project, "script/scan_tissue_attempt/get_gene_annotations.R")))
  writeLines(c("A", "B", "C"), file.path(project, "data/genes_protein_coding.txt"))
  bfile <- file.path(datadir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv")
  for (path in c(paste0(bfile, c(".bed", ".bim", ".fam")), file.path(datadir, "plink2")))
    writeLines("fixture", path)
  Sys.chmod(file.path(datadir, "plink2"), "0755")
  connection <- gzfile(file.path(datadir,
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"), "wt")
  writeLines(c("# fixture",
    'chr1\tensembl\tgene\t1000000\t1010000\t.\t+\t.\tgene_id "ENSA"; gene_type "protein_coding"; gene_name "A";',
    'chr2\tensembl\tgene\t2000000\t2010000\t.\t-\t.\tgene_id "ENSB"; gene_type "protein_coding"; gene_name "B";',
    'chrX\tensembl\tgene\t3000000\t3010000\t.\t+\t.\tgene_id "ENSC"; gene_type "protein_coding"; gene_name "C";'), connection)
  close(connection)

  variables <- c("SUSIE_MIX_GTEX_DIR", "SUSIE_MIX_PLINK", "SUSIE_MIX_GENOTYPE_DIR",
                 "SUSIE_MIX_SIM_GENE", "SUSIE_MIX_PLINK_MEMORY_MB")
  previous <- Sys.getenv(variables, unset = NA_character_)
  on.exit({
    Sys.unsetenv(variables)
    values <- previous[!is.na(previous)]
    if (length(values)) do.call(Sys.setenv, as.list(values))
  }, add = TRUE)
  Sys.unsetenv(variables)
  Sys.setenv(SUSIE_MIX_GTEX_DIR = datadir, SUSIE_MIX_PLINK_MEMORY_MB = "8000")
  empty <- file.path(project, "temp_plink") # deliberately absent
  genotype_source <- prepare_simulation_genotypes(project, empty)
  stopifnot(genotype_source$mode == "gtex", identical(genotype_source$loci$gene, c("A", "B")),
            genotype_source$plink_memory == 8000L, !dir.exists(empty))

  set.seed(9014)
  X <- sapply(seq(.15, .45, length.out = 16), function(maf) rbinom(160, 2, maf))
  colnames(X) <- paste0("snp", seq_len(ncol(X)))
  raw <- data.frame(FID = 1:160, IID = paste0("donor", 1:160), PAT = 0, MAT = 0,
                    SEX = 0, PHENOTYPE = -9, X)
  raw_file <- file.path(raw_dir, "fixture.raw")
  write.table(raw, raw_file, row.names = FALSE, quote = FALSE, sep = "\t")
  cached_source <- prepare_simulation_genotypes(project, raw_dir)
  stopifnot(cached_source$mode == "raw")

  # Missing input files are diagnosed before submission/fitting.
  Sys.setenv(SUSIE_MIX_GTEX_DIR = file.path(directory, "missing"))
  missing <- tryCatch(prepare_simulation_genotypes(project, empty), error = conditionMessage)
  stopifnot(grepl("GTEx extraction needs", missing, fixed = TRUE),
            prepare_simulation_genotypes(project, raw_dir)$mode == "raw")
  Sys.setenv(SUSIE_MIX_GTEX_DIR = datadir)

  produced <- character()
  attempts <- 0L
  fail_on <- 2L
  original_call <- .init_plink_call
  on.exit(assign(".init_plink_call", original_call, envir = .GlobalEnv), add = TRUE)
  assign(".init_plink_call", function(executable, args, logfile) {
    attempts <<- attempts + 1L
    prefix <- args[match("--out", args) + 1L]
    produced <<- c(produced, paste0(prefix, ".raw"))
    stopifnot(file.copy(raw_file, paste0(prefix, ".raw")))
    if (attempts == fail_on) {
      writeLines("Error: Out of memory.", logfile)
      return(2L)
    }
    writeLines("mock extraction succeeded", logfile)
    0L
  }, envir = .GlobalEnv)

  job <- read.csv("script/sim/jobs_slide/manifest.csv")[1L, , drop = FALSE]
  job$n <- 120L
  job$L <- 2L
  job$reps_per_chunk <- 2L
  job$output_file <- "results/test.RData"
  run <- function() run_simulation_job(job$job_id, project, empty, job = job)
  failure <- tryCatch(run(), error = identity)
  stopifnot(inherits(failure, "init_plink_memory_error"), attempts == 2L,
            !any(file.exists(produced)))
  checkpoint <- file.path(project, job$output_file)
  saved <- new.env()
  load(checkpoint, saved)
  first <- saved$results[[1L]]
  stopifnot(length(saved$results) == 1L, is.null(first$error),
            identical(first$metrics$method, sim_methods),
            first$genotype_mode == "gtex_fallback", first$gene %in% c("A", "B"),
            saved$checkpoint_settings$genotype_inputs$mode == "gtex_fallback",
            is.null(saved$checkpoint_settings$arguments$temp_dir))

  # Increasing PLINK memory resumes exactly the unfinished seed, preserving
  # the first real five-method fit and its selected gene.
  Sys.setenv(SUSIE_MIX_PLINK_MEMORY_MB = "16000")
  run()
  load(checkpoint, saved)
  stopifnot(attempts == 3L, length(saved$results) == 2L,
            identical(first, saved$results[[1L]]),
            saved$checkpoint_settings$genotype_inputs$plink_memory == 16000L,
            !any(file.exists(produced)))
  second <- saved$results[[2L]]
  stopifnot(is.null(second$error), identical(second$metrics$method, sim_methods))

  # Reconstruct the same phenotype/causal SNPs through the original cached
  # path and compare all five models, independently of the extraction wrapper.
  arguments <- saved$checkpoint_settings$arguments
  arguments$seed <- second$seed
  arguments$temp_dir <- raw_dir
  direct <- do.call(sim_mix, arguments)
  for (field in c("metrics", "causal_snps", "true_pos", "beta_standardized", sim_pip_fields))
    stopifnot(isTRUE(all.equal(direct[[field]], second[[field]], tolerance = 1e-12)))
  set.seed(second$seed)
  expected_gene <- genotype_source$loci$gene[sample.int(nrow(genotype_source$loci), 1L)]
  stopifnot(identical(expected_gene, second$gene))

  # PVE and generating coding changes preserve the seeded gene choice.
  arguments$return_data <- TRUE
  arguments$pve <- .10
  a <- simulate_slide_replicate(arguments, genotype_source)
  arguments$L_pdom <- 0L
  arguments$L_add <- 1L
  b <- simulate_slide_replicate(arguments, genotype_source)
  stopifnot(identical(a$gene, second$gene), identical(b$gene, second$gene),
            identical(a$X, b$X), !any(file.exists(produced)))
  arguments$n <- 1000L
  bad <- simulate_slide_replicate(arguments, genotype_source)
  stopifnot(grepl("Not enough donors", bad$error), !any(file.exists(produced)))

  before <- tools::md5sum(checkpoint)
  attempts_before <- attempts
  run()
  stopifnot(identical(before, tools::md5sum(checkpoint)), attempts == attempts_before)
  Sys.setenv(SUSIE_MIX_SIM_GENE = "B")
  incompatible <- tryCatch(run(), error = conditionMessage)
  stopifnot(grepl("Checkpoint design", incompatible), attempts == attempts_before)
  Sys.unsetenv("SUSIE_MIX_SIM_GENE")
  changed <- tryCatch(run_simulation_job(job$job_id, project, raw_dir, job = job),
                      error = conditionMessage)
  stopifnot(grepl("Checkpoint design", changed))
  writeLines("new region", file.path(raw_dir, "new.raw"))
  changed <- tryCatch(simulate_slide_replicate(arguments, cached_source), error = conditionMessage)
  stopifnot(grepl("file list changed", changed, fixed = TRUE), file.exists(raw_file))
  cat("PASS: full runner GTEx fallback, five real fits, seeded genes, cleanup, PLINK failure, resume, source validation and cached-pool precedence.\n")
}

test_slide_genotypes()
