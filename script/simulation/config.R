# Source from run_simulation.R, or override entries in an R file containing `config`.
simulation_config <- function(root) {
  datadir <- Sys.getenv("GTEX_DATADIR", "/project2/mstephens/gtex")
  list(
    root = normalizePath(root, winslash = "/", mustWork = TRUE),
    output_dir = file.path(root, "simulation results", "production"),
    genotype_prefix = file.path(datadir,
      "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"),
    plink_exec = Sys.getenv("PLINK2", file.path(datadir, "plink2")),
    gtf_file = file.path(datadir,
      "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"),
    gene_list = file.path(root, "data", "genes_protein_coding.txt"),
    n = 200L, replications = 200L, seed = 20260908L,
    pve = c(0.05, 0.10, 0.20),
    # These genotype and fitting settings match workhorse.R.
    cis_window = 5e5, min_maf_plink = 0, min_maf = 0.05,
    hwe_thresh = 1e-8, min_n_rec = 5L,
    L = 10L, standardize = FALSE, estimate_prior_method = "EM",
    min_abs_corr = 0, coverage = 0.95,
    # Allow more iterations than the workhorse's implicit susieR default.
    max_iter = 1000L, tol = 1e-3, prior_tol = 1e-9, n_purity = 100L,
    pip_threshold = 0.95,
    thresholds = sort(unique(c(seq(0, 1, by = 0.01), 0.95, 0.99, 0.999))),
    calibration_breaks = seq(0, 1, by = 0.05),
    # Normal residuals are independent of genotypes. This targets conditional PVE;
    # realized sample PVE fluctuates and is recorded, not forced by orthogonalization.
    effect_distribution = "equal_standardized_random_sign",
    max_locus_attempts = 100L, max_causal_attempts = 1000L,
    plink_threads = 2L, plink_memory_mb = 8000L,
    save_fits = FALSE, pilot = FALSE, raw_file = NULL
  )
}
