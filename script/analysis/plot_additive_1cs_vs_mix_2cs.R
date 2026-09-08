# plot_additive_1cs_vs_mix_2cs.R
#
# For every significant, sufficiently expressed gene-tissue pair for which
# additive SuSiE has one credible set and SuSiE-mix has two credible sets,
# create a four-panel figure:
#
#   1. additive SuSiE PIPs;
#   2. SuSiE-mix PIPs, with coding-block boundaries;
#   3. expression versus the CS1 lead SNP, adjusted for the CS2 lead predictor;
#   4. expression versus the CS2 lead SNP, adjusted for the CS1 lead predictor.
#
# The phenotype, genotype QC, mixed codings, and tissue-specific predictor
# filtering are reconstructed with the same rules as the current workhorse.

library(data.table)
library(matrixStats)
library(susieR)


# ============================================================
# User settings
# ============================================================

project_dir <- "/project2/mstephens/wdenault/susie_mix"
datadir <- "/project2/mstephens/gtex"

summary_file <- file.path(project_dir, "res_summary.RData")
results_dir <- file.path(project_dir, "results")
temp_dir <- file.path(project_dir, "temp_plink")

plot_dir <- file.path(
  project_dir,
  "plot",
  "additive_1cs_vs_mix_2cs"
)

plot_summary_file <- file.path(
  plot_dir,
  "additive_1cs_vs_mix_2cs_plot_summary.csv"
)

association_threshold <- 1e-8
minimum_mean_reads <- 100

# These values match scan_tissue_attempt/workhorse.R.
min_maf_plink <- 0.00
min_maf <- 0.05
hwe_thresh <- 1e-8
min_n_rec <- 5
cis_window <- 5e5
min_samples <- 50


# ============================================================
# Shared workhorse-compatible reconstruction functions
# ============================================================

source(
  file.path(
    project_dir,
    "script",
    "analysis",
    "fit_mix_expression_plot_utils.R"
  )
)


# ============================================================
# Plotting helpers
# ============================================================

add_mix_coding_boundaries <- function(predictor_map) {

  required_columns <- c("predictor_index", "coding")
  missing_columns <- setdiff(required_columns, names(predictor_map))

  if (length(missing_columns) > 0L) {
    stop(
      "mix_predictor_map is missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  ordered_map <- predictor_map[
    order(predictor_map$predictor_index),
    ,
    drop = FALSE
  ]

  coding <- as.character(ordered_map$coding)

  if (length(coding) == 0L) {
    return(invisible(NULL))
  }

  coding_runs <- rle(coding)
  run_ends <- cumsum(coding_runs$lengths)
  run_starts <- c(1L, head(run_ends, -1L) + 1L)
  run_midpoints <- (run_starts + run_ends) / 2

  if (length(run_ends) > 1L) {
    abline(
      v = head(run_ends, -1L) + 0.5,
      col = "#B22222",
      lty = 2,
      lwd = 1.4
    )
  }

  plot_limits <- par("usr")

  text(
    x = run_midpoints,
    y = plot_limits[4],
    labels = tools::toTitleCase(coding_runs$values),
    pos = 3,
    offset = 0.20,
    xpd = NA,
    cex = 0.68,
    col = "#7A1F1F"
  )

  invisible(NULL)
}


finite_plot_limits <- function(...) {

  values <- unlist(list(...), use.names = FALSE)
  values <- values[is.finite(values)]

  if (length(values) == 0L) {
    stop("The conditional expression values are all non-finite.")
  }

  limits <- range(values)
  spread <- diff(limits)
  padding <- if (spread > 0) 0.06 * spread else 0.25

  limits + c(-padding, padding)
}


open_four_panel_png <- function(filename) {

  png_arguments <- list(
    filename = filename,
    width = 3200,
    height = 2800,
    res = 250
  )

  if (capabilities("cairo")) {
    png_arguments$type <- "cairo-png"
  }

  do.call(png, png_arguments)
}


plot_four_panel_case <- function(
    gene_name,
    tissue_name,
    fit_add,
    fit_mix,
    predictor_map,
    leads,
    raw_genotype_cs1,
    raw_genotype_cs2,
    conditional_expression_cs1,
    conditional_expression_cs2,
    output_file) {

  common_expression_limits <- finite_plot_limits(
    conditional_expression_cs1,
    conditional_expression_cs2
  )

  open_four_panel_png(output_file)
  on.exit(invisible(dev.off()), add = TRUE)

  par(
    mfrow = c(2, 2),
    oma = c(1, 1, 4, 1),
    mar = c(5, 5, 4.5, 2) + 0.1,
    las = 1
  )

  susie_plot(
    fit_add,
    y = "PIP",
    main = "Additive SuSiE: 1 credible set"
  )

  susie_plot(
    fit_mix,
    y = "PIP",
    main = "SuSiE-mix: 2 credible sets"
  )

  add_mix_coding_boundaries(predictor_map)

  cs1_title <- paste0(
    leads$cs_name[1],
    " lead: ",
    leads$lead_snp[1],
    " (",
    leads$lead_coding[1],
    ")\nAdjusted for ",
    leads$cs_name[2],
    " lead (",
    leads$lead_coding[2],
    ")"
  )

  plot_fit_mix_expression_panel(
    lead = leads[1],
    raw_genotype = raw_genotype_cs1,
    expression_to_plot = conditional_expression_cs1,
    main = cs1_title,
    ylab = "Conditional normalized gene expression",
    ylim = common_expression_limits,
    show_mean_legend = FALSE,
    jitter_seed = 101L
  )

  cs2_title <- paste0(
    leads$cs_name[2],
    " lead: ",
    leads$lead_snp[2],
    " (",
    leads$lead_coding[2],
    ")\nAdjusted for ",
    leads$cs_name[1],
    " lead (",
    leads$lead_coding[1],
    ")"
  )

  plot_fit_mix_expression_panel(
    lead = leads[2],
    raw_genotype = raw_genotype_cs2,
    expression_to_plot = conditional_expression_cs2,
    main = cs2_title,
    ylab = "Conditional normalized gene expression",
    ylim = common_expression_limits,
    show_mean_legend = TRUE,
    jitter_seed = 102L
  )

  mtext(
    paste(gene_name, tissue_name, sep = " - "),
    side = 3,
    outer = TRUE,
    line = 1.4,
    font = 2,
    cex = 1.25
  )

  invisible(NULL)
}


# ============================================================
# Select the same analysis population as res_idx
# ============================================================

summary_environment <- new.env(parent = emptyenv())
loaded_objects <- load(summary_file, envir = summary_environment)

if ("res_summary" %in% loaded_objects) {
  res <- as.data.table(summary_environment$res_summary)
} else {
  stop("res_summary.RData does not contain an object named res_summary.")
}

required_summary_columns <- c(
  "gene",
  "tissue",
  "min_pv",
  "mean_count",
  "ncs_susie",
  "ncs_susie_mix",
  "overlap_snp"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  names(res)
)

if (length(missing_summary_columns) > 0L) {
  stop(
    "res_summary is missing: ",
    paste(missing_summary_columns, collapse = ", "),
    "."
  )
}

cases <- unique(
  res[
    is.finite(min_pv) &
      min_pv < association_threshold &
      is.finite(mean_count) &
      mean_count >= minimum_mean_reads &
      !is.na(ncs_susie) &
      !is.na(ncs_susie_mix) &
      !is.na(overlap_snp) &
      ncs_susie == 1L &
      ncs_susie_mix == 2L,
    .(
      gene,
      tissue,
      min_pv,
      mean_count,
      ncs_susie,
      ncs_susie_mix
    )
  ],
  by = c("gene", "tissue")
)

setorder(cases, gene, tissue)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(cases) == 0L) {
  stop("No additive-1-CS / SuSiE-mix-2-CS cases met the filters.")
}

cat(
  "Selected ",
  nrow(cases),
  " gene-tissue pairs (",
  uniqueN(cases$gene),
  " genes).\n",
  sep = ""
)


# ============================================================
# Load shared GTEx inputs once
# ============================================================

plink_exec <- file.path(datadir, "plink2")
geno_file <- file.path(
  datadir,
  "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"
)
subject_pheno_file <- file.path(
  datadir,
  "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"
)
sample_attr_file <- file.path(
  datadir,
  "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"
)
expr_file <- file.path(
  datadir,
  "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
)
gtf_file <- file.path(
  datadir,
  "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"
)
gene_annotation_function <- file.path(
  project_dir,
  "script",
  "scan_tissue_attempt",
  "get_gene_annotations.R"
)

shared_inputs <- load_fit_mix_expression_inputs(
  subject_pheno_file = subject_pheno_file,
  sample_attr_file = sample_attr_file,
  expr_file = expr_file,
  gtf_file = gtf_file,
  gene_annotation_function = gene_annotation_function
)


# ============================================================
# Reconstruct each gene once, then plot all selected tissues
# ============================================================

result_rows <- list()
result_index <- 0L

add_result <- function(result) {
  result_index <<- result_index + 1L
  result_rows[[result_index]] <<- result

  # Keep an up-to-date audit file even if an RCC allocation ends before
  # every case has been processed.
  fwrite(
    rbindlist(result_rows, fill = TRUE),
    plot_summary_file
  )
}

make_error_row <- function(gene_name, tissue_name, message) {
  data.table(
    gene = gene_name,
    tissue = tissue_name,
    status = "error",
    message = message
  )
}

for (gene_name in unique(cases$gene)) {

  gene_cases <- cases[gene == gene_name]

  gene_data <- tryCatch(
    prepare_fit_mix_gene_data(
      gene_name = gene_name,
      shared_inputs = shared_inputs,
      results_dir = results_dir,
      temp_dir = temp_dir,
      plink_exec = plink_exec,
      geno_file = geno_file,
      cis_window = cis_window,
      min_maf_plink = min_maf_plink,
      min_maf = min_maf,
      hwe_thresh = hwe_thresh
    ),
    error = function(e) e
  )

  if (inherits(gene_data, "error")) {
    error_message <- conditionMessage(gene_data)

    for (tissue_name in gene_cases$tissue) {
      message(
        "ERROR [",
        gene_name,
        " / ",
        tissue_name,
        "]: ",
        error_message
      )
      add_result(make_error_row(gene_name, tissue_name, error_message))
    }

    next
  }

  for (tissue_name in gene_cases$tissue) {

    case_result <- tryCatch(
      {
        if (!tissue_name %in% names(gene_data$gene_results)) {
          stop("Tissue is absent from the saved gene result.")
        }

        tissue_result <- gene_data$gene_results[[tissue_name]]
        fit_add <- tissue_result$susie_add
        fit_mix <- tissue_result$susie_mix

        if (is.null(fit_add)) {
          stop("The saved tissue result has no susie_add fit.")
        }

        if (is.null(fit_mix)) {
          stop("The saved tissue result has no susie_mix fit.")
        }

        add_cs <- fit_add$sets$cs
        mix_cs <- fit_mix$sets$cs
        n_add_cs <- if (is.null(add_cs)) 0L else length(add_cs)
        n_mix_cs <- if (is.null(mix_cs)) 0L else length(mix_cs)

        if (n_add_cs != 1L || n_mix_cs != 2L) {
          stop(
            "Summary expected 1 additive CS and 2 mixed CSs, but the saved fits have ",
            n_add_cs,
            " and ",
            n_mix_cs,
            "."
          )
        }

        tissue_data <- prepare_fit_mix_tissue_data(
          gene_data = gene_data,
          gene_name = gene_name,
          tissue_name = tissue_name,
          tissue_result = tissue_result,
          min_samples = min_samples,
          min_n_rec = min_n_rec
        )

        leads <- get_fit_mix_cs_leads(
          fit_mix = fit_mix,
          predictor_map = tissue_data$predictor_map
        )

        if (nrow(leads) != 2L) {
          stop("Could not identify exactly two SuSiE-mix CS leads.")
        }

        missing_lead_snps <- setdiff(
          leads$lead_snp,
          colnames(tissue_data$geno_for_counts)
        )

        if (length(missing_lead_snps) > 0L) {
          stop(
            "Lead SNPs are absent from reconstructed genotypes: ",
            paste(missing_lead_snps, collapse = ", "),
            "."
          )
        }

        raw_genotype_cs1 <- tissue_data$geno_for_counts[
          ,
          leads$lead_snp[1]
        ]
        raw_genotype_cs2 <- tissue_data$geno_for_counts[
          ,
          leads$lead_snp[2]
        ]

        conditional_expression_cs1 <- conditional_expression_for_cs(
          y = tissue_data$y,
          geno_mix = tissue_data$geno_mix,
          lead_indices = leads$lead_predictor_index,
          focal_cs = 1L
        )

        conditional_expression_cs2 <- conditional_expression_for_cs(
          y = tissue_data$y,
          geno_mix = tissue_data$geno_mix,
          lead_indices = leads$lead_predictor_index,
          focal_cs = 2L
        )

        output_file <- file.path(
          plot_dir,
          paste0(
            sanitize_filename(gene_name),
            "_",
            sanitize_filename(tissue_name),
            "_add1_mix2.png"
          )
        )

        plot_four_panel_case(
          gene_name = gene_name,
          tissue_name = tissue_name,
          fit_add = fit_add,
          fit_mix = fit_mix,
          predictor_map = tissue_data$predictor_map,
          leads = leads,
          raw_genotype_cs1 = raw_genotype_cs1,
          raw_genotype_cs2 = raw_genotype_cs2,
          conditional_expression_cs1 = conditional_expression_cs1,
          conditional_expression_cs2 = conditional_expression_cs2,
          output_file = output_file
        )

        data.table(
          gene = gene_name,
          tissue = tissue_name,
          status = "plotted",
          message = "Four-panel plot created",
          additive_n_cs = n_add_cs,
          mix_n_cs = n_mix_cs,
          cs1_name = leads$cs_name[1],
          cs1_component = leads$component_index[1],
          cs1_lead_predictor = leads$lead_predictor_name[1],
          cs1_lead_snp = leads$lead_snp[1],
          cs1_lead_coding = leads$lead_coding[1],
          cs1_lead_pip = leads$lead_pip[1],
          cs1_size = leads$cs_size[1],
          cs1_n_genotype_0 = sum(raw_genotype_cs1 == 0L),
          cs1_n_genotype_1 = sum(raw_genotype_cs1 == 1L),
          cs1_n_genotype_2 = sum(raw_genotype_cs1 == 2L),
          cs1_adjusted_for = leads$lead_predictor_name[2],
          cs2_name = leads$cs_name[2],
          cs2_component = leads$component_index[2],
          cs2_lead_predictor = leads$lead_predictor_name[2],
          cs2_lead_snp = leads$lead_snp[2],
          cs2_lead_coding = leads$lead_coding[2],
          cs2_lead_pip = leads$lead_pip[2],
          cs2_size = leads$cs_size[2],
          cs2_n_genotype_0 = sum(raw_genotype_cs2 == 0L),
          cs2_n_genotype_1 = sum(raw_genotype_cs2 == 1L),
          cs2_n_genotype_2 = sum(raw_genotype_cs2 == 2L),
          cs2_adjusted_for = leads$lead_predictor_name[1],
          output_file = output_file
        )
      },
      error = function(e) {
        error_message <- conditionMessage(e)
        message(
          "ERROR [",
          gene_name,
          " / ",
          tissue_name,
          "]: ",
          error_message
        )
        make_error_row(gene_name, tissue_name, error_message)
      }
    )

    add_result(case_result)
  }
}

plot_summary <- rbindlist(result_rows, fill = TRUE)
fwrite(plot_summary, plot_summary_file)

cat(
  "\nFinished four-panel plots. Plotted ",
  sum(plot_summary$status == "plotted"),
  " of ",
  nrow(plot_summary),
  " selected gene-tissue pairs.\n",
  "Plots: ",
  plot_dir,
  "\nSummary: ",
  plot_summary_file,
  "\n",
  sep = ""
)
