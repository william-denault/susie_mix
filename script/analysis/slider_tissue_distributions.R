# Compare distributions of reported posterior means, conditional on selection.
# Source this file, then call slider_tissue_distributions(results).
# No packages are required; 'cluster' enables silhouette summaries if installed.
#
# Each file must identify one gene. All records for a resampling block stay
# together, across tissues and CSs. Defaults: equal weight per gene within each
# tissue, unforced leads only, and all available genes in each tissue.
# These are different estimands from a distribution of true biological deltas.
#
# Inference: an asymptotic, centered cluster bootstrap of CDF contrasts. It does
# not permute tissue labels. Blocks must be approximately independent and
# representative of the target gene population. Use a locus/block column when
# dependence across nearby genes is material. Fitted posteriors and selection
# are treated as the observed outcomes; priors are not re-estimated here.
#
# CDFs are evaluated on a fixed grid. W1 distances and integrated squared-CDF
# statistics are numerical approximations; increase n_grid to check resolution.
# Pairwise tests use each tissue's selected cohort, NOT pairwise common genes.
# For a matched-cohort comparison subset to the tissues of interest and use
# cohort = "common". This still conditions on detection in those tissues.
#
# Example:
# a <- slider_tissue_distributions(results, B = 9999, seed = 1)
# a$omnibus
# a$pairwise[order(a$pairwise$p_holm), ]
# a$sample_sizes
# a$silhouette
# plot(a$tree, xlab = "", sub = "Wasserstein-1; average linkage")
# a$clusters
# heatmap(a$coassignment[["3"]], symm = TRUE, scale = "none")
#
# Cluster stability is descriptive: a high coassignment rate does not establish
# that there are discrete biological classes. Compare k = 2,...,6, including
# the user's candidate k = 3 and 4. Never select k by the omnibus p-value.

slider_tissue_distributions <- function(
    results, value_col = "lead_delta_mean", gene_col = "file",
    tissue_col = "tissue", forced_col = "lead_delta_forced",
    block_col = gene_col, weighting = c("gene", "signal"),
    cohort = c("available", "common"), B = 9999L, n_grid = 501L,
    ks = 2:6, seed = 1L, verbose = TRUE) {
  weighting <- match.arg(weighting)
  cohort <- match.arg(cohort)
  required <- unique(c(value_col, gene_col, tissue_col, block_col, forced_col))
  if (!is.data.frame(results) || !all(required %in% names(results)))
    stop("results must contain: ", paste(required, collapse = ", "))
  if (!is.numeric(results[[value_col]])) stop("value_col must be numeric.")
  if (length(B) != 1L || !is.finite(B) || B < 99 || B != as.integer(B))
    stop("B must be an integer >= 99.")
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 21 ||
      n_grid != as.integer(n_grid)) stop("n_grid must be an integer >= 21.")
  forced <- if (is.null(forced_col)) rep(FALSE, nrow(results)) else
    results[[forced_col]]
  if (!is.logical(forced)) stop("forced_col must be logical (or NULL).")
  d <- data.frame(x = results[[value_col]],
                  gene = as.character(results[[gene_col]]),
                  tissue = as.character(results[[tissue_col]]),
                  block = as.character(results[[block_col]]))
  keep <- !is.na(forced) & !forced & is.finite(d$x) &
    !is.na(d$gene) & nzchar(d$gene) & !is.na(d$tissue) & nzchar(d$tissue) &
    !is.na(d$block) & nzchar(d$block)
  n_filtered <- sum(!keep)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("No valid unforced observations.")
  if (any(abs(d$x) > 1 + 1e-10)) stop("Posterior means must be in [-1, 1].")
  d$x <- pmax(-1, pmin(1, d$x))
  tissues <- sort(unique(d$tissue))
  nt <- length(tissues)
  if (nt < 2L) stop("Need at least two tissues.")
  if (choose(nt, 2L) / (B + 1) > 0.05)
    warning("B is too small for any Holm-adjusted pairwise p-value below 0.05; ",
            "increase B for pairwise inference. The omnibus test is still usable.")
  if (cohort == "common") {
    present <- unique(d[c("gene", "tissue")])
    counts <- table(present$gene)
    common <- names(counts)[counts == nt]
    d <- d[d$gene %in% common, , drop = FALSE]
    if (!nrow(d)) stop("No eligible genes shared by all selected tissues.")
  }
  gene_blocks <- unique(d[c("gene", "block")])
  if (anyDuplicated(gene_blocks$gene))
    stop("Every gene must belong to the same block across all tissues.")
  genes <- unique(d$gene)
  blocks <- unique(d$block)
  nb <- length(blocks)
  if (nb < 2L) stop("Need at least two independent resampling blocks.")
  ti <- match(d$tissue, tissues)
  gi <- match(d$gene, genes)
  bi <- match(d$block, blocks)
  gene_tissue <- (gi - 1L) * nt + ti
  weights <- if (weighting == "gene")
    1 / tabulate(gene_tissue)[gene_tissue] else rep(1, nrow(d))
  grid <- seq(-1, 1, length.out = n_grid)
  dx <- grid[2L] - grid[1L]
  # First grid knot >= x. Thus the CDF at each knot includes x equal to it.
  bin <- findInterval(d$x, grid, left.open = TRUE) + 1L
  bin <- pmax(1L, pmin(n_grid, bin))
  cell <- (ti - 1L) * n_grid + bin
  occupied <- sort(unique(cell))
  make_cdf <- function(row_weights) {
    sums <- rowsum(row_weights, cell, reorder = TRUE)[, 1L]
    hist <- numeric(nt * n_grid)
    hist[occupied] <- sums
    hist <- matrix(hist, nrow = nt, byrow = TRUE)
    cumulative <- t(apply(hist, 1L, cumsum))
    mass <- cumulative[, n_grid]
    if (any(mass <= 0)) return(NULL)
    sweep(cumulative, 1L, mass, "/")
  }
  F <- make_cdf(weights)
  rownames(F) <- tissues
  # Exclude the last knot: the CDF is identically one there.
  integral_cols <- seq_len(n_grid - 1L)
  global_stat <- function(curves) {
    centered <- sweep(curves[, integral_cols, drop = FALSE], 2L,
                      colMeans(curves[, integral_cols, drop = FALSE]), "-")
    sum(centered^2) * dx
  }
  pair_stat <- function(curves)
    as.numeric(dist(curves[, integral_cols, drop = FALSE]))^2 * dx
  w1_dist <- function(curves)
    dist(curves[, integral_cols, drop = FALSE], method = "manhattan") * dx
  observed_global <- global_stat(F)
  observed_pairs <- pair_stat(F)
  W1 <- w1_dist(F)
  tree <- hclust(W1, method = "average")
  ks <- sort(unique(as.integer(ks)))
  ks <- ks[is.finite(ks) & ks >= 2L & ks < nt]
  clusters <- lapply(ks, function(k) cutree(tree, k = k))
  names(clusters) <- as.character(ks)
  coassignment <- lapply(ks, function(k)
    matrix(0, nt, nt, dimnames = list(tissues, tissues)))
  names(coassignment) <- names(clusters)
  sample_sizes <- do.call(rbind, lapply(tissues, function(tissue) {
    dt <- d[d$tissue == tissue, , drop = FALSE]
    data.frame(tissue = tissue, n_genes = length(unique(dt$gene)),
               n_blocks = length(unique(dt$block)), n_signals = nrow(dt))
  }))
  if (any(sample_sizes$n_blocks < 30L))
    warning("Some tissues have fewer than 30 blocks; bootstrap inference is fragile.")
  set.seed(seed)
  exceed_global <- 0L
  exceed_pairs <- integer(length(observed_pairs))
  used <- 0L
  attempts <- 0L
  while (used < B && attempts < 2L * B) {
    attempts <- attempts + 1L
    multiplicity <- tabulate(sample.int(nb, nb, replace = TRUE), nbins = nb)
    Fb <- make_cdf(weights * multiplicity[bi])
    if (is.null(Fb)) next
    used <- used + 1L
    # Center the CDF contrasts before obtaining null critical values.
    # Comparing uncentered bootstrap distances with observed distances is NOT
    # a test of equality: ordinary resampling retains the observed difference.
    null_curves <- Fb - F
    exceed_global <- exceed_global + (global_stat(null_curves) >= observed_global)
    exceed_pairs <- exceed_pairs + (pair_stat(null_curves) >= observed_pairs)
    if (length(ks)) {
      boot_tree <- hclust(w1_dist(Fb), method = "average")
      for (j in seq_along(ks)) {
        labels <- cutree(boot_tree, k = ks[j])
        coassignment[[j]] <- coassignment[[j]] + outer(labels, labels, "==")
      }
    }
    if (verbose && used %% 250L == 0L) message(used, "/", B, " bootstrap draws")
  }
  if (used < B) stop("Too many bootstrap draws with empty tissues; use larger cohorts.")
  if (attempts > B)
    warning("Skipped ", attempts - B, " empty-tissue draws; sparse cohorts need care.")
  coassignment <- lapply(coassignment, function(x) x / B)
  pair_index <- which(lower.tri(matrix(0, nt, nt)), arr.ind = TRUE)
  pairwise <- data.frame(tissue1 = tissues[pair_index[, 2L]],
                         tissue2 = tissues[pair_index[, 1L]],
                         W1 = as.numeric(W1),
                         integrated_squared_cdf = observed_pairs,
                         p_value = (1 + exceed_pairs) / (B + 1))
  # Holm controls familywise error without independence of the pairwise tests.
  pairwise$p_holm <- p.adjust(pairwise$p_value, method = "holm")
  silhouette <- data.frame(k = ks, mean_silhouette = rep(NA_real_, length(ks)))
  if (length(ks) && requireNamespace("cluster", quietly = TRUE) && any(W1 > 0))
    silhouette$mean_silhouette <- vapply(clusters, function(labels)
      mean(cluster::silhouette(labels, W1)[, "sil_width"]), numeric(1L))
  list(omnibus = data.frame(statistic = observed_global,
                            p_value = (1 + exceed_global) / (B + 1), B = B),
       pairwise = pairwise, distance = W1, tree = tree,
       clusters = clusters, silhouette = silhouette,
       coassignment = coassignment, sample_sizes = sample_sizes,
       cdf = F, grid = grid,
       settings = list(weighting = weighting, cohort = cohort,
                       value_col = value_col, block_col = block_col,
                       n_filtered = n_filtered, n_blocks = nb,
                       n_grid = n_grid, seed = seed,
                       minimum_p = 1 / (B + 1)))
}
