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
# Each title also reports the one-CS plots' ELBO + KL likelihood comparison
# and min(|additive lead - mixed CS1 lead|, |additive lead - mixed CS2 lead|)
# in GRCh38 base pairs. The audit CSV retains all three pairwise distances.

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

get_add1_mix2_lead_distances <- function(add_lead, leads) {

  if (length(add_lead$lead_snp) != 1L || nrow(leads) != 2L) {
    stop("Expected one additive lead and two mixed CS leads.")
  }
  ids <- as.character(c(add_lead$lead_snp, leads$lead_snp))
  valid <- !is.na(ids) & grepl(
    "^chr[^_]+_[0-9]+_[^_]+_[^_]+_b38(?:_|$)", ids, perl = TRUE
  )
  chromosomes <- rep(NA_character_, 3L)
  positions <- rep(NA_real_, 3L)
  chromosomes[valid] <- sub("^(chr[^_]+)_.*$", "\\1", ids[valid])
  positions[valid] <- as.numeric(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", ids[valid]))
  valid <- valid & is.finite(positions) & positions > 0
  distance <- function(i, j) {
    if (valid[i] && valid[j] && chromosomes[i] == chromosomes[j]) {
      abs(positions[i] - positions[j])
    } else {
      NA_real_
    }
  }
  additive_distances <- c(distance(1L, 2L), distance(1L, 3L))
  # Both distances are needed to establish the minimum. Do not silently
  # choose a winner when one coordinate pair is unavailable.
  comparable <- all(is.finite(additive_distances))
  closest_distance <- if (comparable) min(additive_distances) else NA_real_
  closest <- if (comparable) which(additive_distances == closest_distance) else integer()

  list(
    additive_lead_chromosome = chromosomes[1],
    additive_lead_position_bp = positions[1],
    cs1_lead_chromosome = chromosomes[2],
    cs1_lead_position_bp = positions[2],
    cs2_lead_chromosome = chromosomes[3],
    cs2_lead_position_bp = positions[3],
    additive_to_cs1_distance_bp = additive_distances[1],
    additive_to_cs2_distance_bp = additive_distances[2],
    cs1_to_cs2_distance_bp = distance(2L, 3L),
    closest_lead_distance_bp = closest_distance,
    closest_mix_cs_names = if (comparable) paste(leads$cs_name[closest], collapse = ";") else NA_character_,
    closest_mix_lead_snps = if (comparable) paste(leads$lead_snp[closest], collapse = ";") else NA_character_,
    closest_mix_lead_codings = if (comparable) paste(leads$lead_coding[closest], collapse = ";") else NA_character_,
    closest_lead_tie = if (comparable) length(closest) > 1L else NA,
    lead_distance_status = if (comparable) "ok" else "unavailable: non-comparable GRCh38 coordinates"
  )
}


plot_four_panel_case <- function(
    gene_name,
    tissue_name,
    fit_add,
    fit_mix,
    predictor_map,
    leads,
    add_lead,
    likelihood_comparison,
    lead_distances,
    raw_genotype_cs1,
    raw_genotype_cs2,
    conditional_expression_cs1,
    conditional_expression_cs2,
    output_file) {

  leads <- as.data.frame(leads)
  statistic <- likelihood_comparison$likelihood_statistic
  statistic_label <- if (is.finite(statistic)) {
    sprintf(
      "Likelihood-ratio statistic: 2(log lik mix - log lik additive) = %.2f",
      statistic
    )
  } else {
    "Likelihood-ratio statistic unavailable (missing finite ELBO/KL metric)"
  }
  distance_label <- if (is.finite(lead_distances$closest_lead_distance_bp)) {
    paste0(
      "Closest lead distance: additive -> ",
      gsub(";", " / ", lead_distances$closest_mix_cs_names, fixed = TRUE),
      " = ", format(lead_distances$closest_lead_distance_bp, big.mark = ",",
                      scientific = FALSE, trim = TRUE),
      " bp (GRCh38)", if (isTRUE(lead_distances$closest_lead_tie)) " [tie]" else ""
    )
  } else {
    "Closest lead distance unavailable (non-comparable GRCh38 coordinates)"
  }

  common_expression_limits <- finite_plot_limits(
    conditional_expression_cs1,
    conditional_expression_cs2
  )

  open_four_panel_png(output_file)
  on.exit(invisible(dev.off()), add = TRUE)

  par(
    mfrow = c(2, 2),
    oma = c(1, 1, 8, 1),
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
    lead = leads[1, , drop = FALSE],
    raw_genotype = raw_genotype_cs1,
    expression_to_plot = conditional_expression_cs1,
    main = cs1_title,
    ylab = "Conditional normalized gene expression",
    ylim = common_expression_limits,
    show_mean_legend = FALSE,
    jitter_seed = 101L,
    genotype_axis_ticks = FALSE
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
    lead = leads[2, , drop = FALSE],
    raw_genotype = raw_genotype_cs2,
    expression_to_plot = conditional_expression_cs2,
    main = cs2_title,
    ylab = "Conditional normalized gene expression",
    ylim = common_expression_limits,
    show_mean_legend = TRUE,
    jitter_seed = 102L,
    genotype_axis_ticks = FALSE
  )

  mtext(
    paste(gene_name, tissue_name, sep = " - "),
    side = 3,
    outer = TRUE,
    line = 6.2,
    font = 2,
    cex = 1.25
  )
  mtext(paste0("Additive lead: ", add_lead$lead_snp),
        side = 3, outer = TRUE, line = 4.8, cex = .85)
  mtext(distance_label, side = 3, outer = TRUE, line = 3.4, cex = .9)
  mtext(statistic_label, side = 3, outer = TRUE, line = 2, cex = .85)
  mtext("ELBO + KL summary metric; positive favors SuSiE-mix; no calibrated p-value",
        side = 3, outer = TRUE, line = .8, cex = .70, col = "gray35")

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

        add_lead <- get_fit_add_plot_lead(
          fit_add = fit_add,
          predictor_map = tissue_result$add_predictor_map
        )
        # This helper uses final ELBO + sum(KL) for either CS count, matching
        # the statistic and unavailable-value handling in the one-CS plots.
        likelihood_comparison <- get_one_cs_likelihood_comparison(fit_add, fit_mix)
        lead_distances <- get_add1_mix2_lead_distances(add_lead, leads)

        missing_lead_snps <- setdiff(
          c(add_lead$lead_snp, leads$lead_snp),
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
          add_lead = add_lead,
          likelihood_comparison = likelihood_comparison,
          lead_distances = lead_distances,
          raw_genotype_cs1 = raw_genotype_cs1,
          raw_genotype_cs2 = raw_genotype_cs2,
          conditional_expression_cs1 = conditional_expression_cs1,
          conditional_expression_cs2 = conditional_expression_cs2,
          output_file = output_file
        )

        case_summary <- data.table(
          gene = gene_name,
          tissue = tissue_name,
          status = "plotted",
          message = "Four-panel plot created",
          additive_n_cs = n_add_cs,
          mix_n_cs = n_mix_cs,
          additive_cs_name = add_lead$cs_name,
          additive_component_index = add_lead$component_index,
          additive_lead_predictor = add_lead$lead_predictor_name,
          additive_lead_snp = add_lead$lead_snp,
          additive_lead_coding = add_lead$lead_coding,
          additive_lead_pip = add_lead$lead_pip,
          additive_cs_size = add_lead$cs_size,
          log_lik_add = likelihood_comparison$log_lik_add,
          log_lik_mix = likelihood_comparison$log_lik_mix,
          likelihood_statistic = likelihood_comparison$likelihood_statistic,
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
        cbind(case_summary, as.data.table(lead_distances))
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
