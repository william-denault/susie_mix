# Select individual CSs whose biological lead SNP is absent from every CS of
# the other model in the SAME gene/tissue. Coding differences alone are not
# TSS disagreement. Shared leads are removed even in partially shared regions.
select_tss_disagreement <- function(cs, additive_model = "SuSiE", mixed_model = "SuSiE-mix",
                                    plot_limit_kb = 200, bin_width_kb = 10) {
  required <- c("gene", "tissue", "model", "lead_snp", "distance_to_tss_kb")
  if (!all(required %in% names(cs))) stop("CS summary is missing TSS agreement fields: ",
                                         paste(setdiff(required, names(cs)), collapse = ", "))
  x <- as.data.frame(cs)
  models <- c(additive_model, mixed_model)
  if (length(models) != 2L || anyNA(models) || anyDuplicated(models)) stop("Specify two distinct models.")
  x <- x[x$model %in% models, , drop = FALSE]
  if (anyNA(x[c("gene", "tissue")]) || any(!nzchar(x$gene) | !nzchar(x$tissue))) {
    stop("TSS comparisons need nonempty gene and tissue names.")
  }
  part <- function(v) paste0(nchar(v), ":", v)
  region <- paste0(part(x$gene), part(x$tissue))
  lead <- trimws(as.character(x$lead_snp))
  valid_lead <- !is.na(lead) & nzchar(lead)
  lead_key <- paste0(region, part(lead))
  x$tss_other_model_cs_count <- integer(nrow(x))
  x$tss_shared_lead <- rep(FALSE, nrow(x))
  x$tss_agreement_status <- rep("unknown_lead", nrow(x))
  for (model_name in models) {
    own <- which(x$model == model_name)
    other <- which(x$model != model_name)
    counts <- table(region[other])
    other_n <- as.integer(counts[region[own]])
    other_n[is.na(other_n)] <- 0L
    shared <- valid_lead[own] & lead_key[own] %in% lead_key[other[valid_lead[other]]]
    unknown_other <- region[own] %in% region[other[!valid_lead[other]]]
    status <- ifelse(shared, "agreed_lead", ifelse(other_n == 0L, "other_model_no_cs",
              ifelse(unknown_other, "unknown_other_model_lead", "different_lead")))
    status[!valid_lead[own]] <- "unknown_lead"
    x$tss_other_model_cs_count[own] <- other_n
    x$tss_shared_lead[own] <- shared
    x$tss_agreement_status[own] <- status
  }
  x$tss_include <- x$tss_agreement_status %in% c("different_lead", "other_model_no_cs")
  # Match the existing centered-bin window, including its outer half bins.
  x$tss_finite_distance <- is.finite(x$distance_to_tss_kb)
  x$tss_in_plot_window <- x$tss_finite_distance &
    abs(x$distance_to_tss_kb) <= plot_limit_kb + bin_width_kb / 2
  x$tss_agreement_rule <- rep("same biological lead SNP within gene/tissue; coding ignored", nrow(x))
  summaries <- list()
  for (tissue_name in c("All tissues", sort(unique(as.character(x$tissue))))) {
    for (model_name in models) {
      z <- x[x$model == model_name & (tissue_name == "All tissues" | x$tissue == tissue_name), , drop = FALSE]
      summaries[[length(summaries) + 1L]] <- data.frame(
        tissue = tissue_name, model = model_name, n_cs_input = nrow(z),
        n_cs_agreed_excluded = sum(z$tss_shared_lead),
        n_cs_different_lead = sum(z$tss_agreement_status == "different_lead"),
        n_cs_other_model_no_cs = sum(z$tss_agreement_status == "other_model_no_cs"),
        n_cs_unknown_excluded = sum(grepl("^unknown", z$tss_agreement_status)),
        n_cs_selected = sum(z$tss_include),
        n_selected_missing_distance = sum(z$tss_include & !z$tss_finite_distance),
        n_selected_outside_window = sum(z$tss_include & z$tss_finite_distance & !z$tss_in_plot_window),
        n_cs_in_window = sum(z$tss_include & z$tss_in_plot_window))
    }
  }
  audit_columns <- unique(c(intersect(c("gene", "tissue", "model", "model_key", "cs_id", "cs_number",
                                         "lead_snp", "lead_coding", "distance_to_tss_kb", "em_iteration"), names(x)),
                            grep("^tss_", names(x), value = TRUE)))
  list(selected = x[x$tss_include, , drop = FALSE], audit = x[audit_columns],
       summary = do.call(rbind, summaries))
}
