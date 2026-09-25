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
  prepare <- function(target_gene = "") prepare_additive_genotypes(
    genotype_dir = file.path(directory, "empty"), datadir = datadir,
    project_dir = project, target_gene = target_gene)
  source <- prepare()
  stopifnot(identical(source$mode, "gtex"), identical(source$loci$gene, c("A", "B")),
    identical(source$loci$chromosome, c("1", "2")), source$plink_threads == 1L,
    identical(as.numeric(source$loci$tss), c(1000000, 2010000)),
    identical(as.numeric(source$loci$from_bp), c(500000, 1510000)),
    identical(prepare("B")$loci$gene, "B"))
  failure <- tryCatch(prepare("Missing"), error = conditionMessage)
  stopifnot(grepl("No matching autosomal", failure, fixed = TRUE))

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
    value("--recode") == "A", value("--maf") == "0", value("--threads") == "1")
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
}
test_gt_extraction()
