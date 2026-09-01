
# susie_cs_overlap.R
#
# Assess overlap between the credible sets (CS) of two susie fit objects.
#
# IMPORTANT: susie CS members are stored as column INDICES into the X
# matrix used for that fit. Comparing indices directly is only valid if
# both fits were run on the same variant ordering (e.g. same gene, same
# genotype matrix, different tissues/phenotypes). If that's not
# guaranteed, pass variant_names1 / variant_names2 (the colnames(X) used
# for each fit) so the function compares by variant ID instead of raw
# index.
#
# Returns a list with:
#   n_cs_1, n_cs_2        -- number of CS found in each fit
#   overlap_matrix        -- n_cs_1 x n_cs_2 matrix, count of shared variants
#   jaccard_matrix         -- same shape, Jaccard index (0-1) per CS pair
#   n_overlapping_cs_1    -- how many of fit1's CS share >=1 variant with any fit2 CS
#   n_overlapping_cs_2    -- same, from fit2's perspective
#   pairs                 -- data.frame of all CS pairs with nonzero overlap

susie_cs_overlap <- function(fit1,fit2,
                             variant_names1 = NULL,
                             variant_names2 = NULL) {

  cs1 <- fit1$sets$cs
  cs2 <- fit2$sets$cs

  n1 <- if (is.null(cs1)) 0 else length(cs1)
  n2 <- if (is.null(cs2)) 0 else length(cs2)

  # Handle the no-CS case(s) cleanly -- e.g. res2$sets$cs is NULL above.
  if (n1 == 0 || n2 == 0) {
    return(list(
      n_cs_1 = n1,n_cs_2 = n2,
      overlap_matrix = NULL,jaccard_matrix = NULL,
      n_overlapping_cs_1 = 0,n_overlapping_cs_2 = 0,
      pairs = data.frame(cs1 = character(0),cs2 = character(0),
                         n_shared_variants = integer(0),jaccard = numeric(0))
    ))
  }

  # Optionally translate indices -> variant names before comparing, so
  # CS from fits on different variant orderings can still be compared
  # correctly (e.g. two tissues where a few SNPs failed QC differently).
  if (!is.null(variant_names1)) cs1 <- lapply(cs1,function(idx) variant_names1[idx])
  if (!is.null(variant_names2)) cs2 <- lapply(cs2,function(idx) variant_names2[idx])

  overlap_n <- matrix(0L,n1,n2,dimnames = list(names(cs1),names(cs2)))
  jaccard   <- matrix(0, n1,n2,dimnames = list(names(cs1),names(cs2)))

  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      inter <- length(intersect(cs1[[i]],cs2[[j]]))
      uni   <- length(union(cs1[[i]],cs2[[j]]))
      overlap_n[i,j] <- inter
      jaccard[i,j]   <- if (uni > 0) inter / uni else 0
    }
  }

  hit_idx  <- which(overlap_n > 0,arr.ind = TRUE)
  pairs_df <- data.frame(
    cs1               = rownames(overlap_n)[hit_idx[,1]],
    cs2               = colnames(overlap_n)[hit_idx[,2]],
    n_shared_variants = overlap_n[hit_idx],
    jaccard            = jaccard[hit_idx]
  )

  list(
    n_cs_1              = n1,
    n_cs_2              = n2,
    overlap_matrix      = overlap_n,
    jaccard_matrix      = jaccard,
    n_overlapping_cs_1  = sum(rowSums(overlap_n > 0) > 0),
    n_overlapping_cs_2  = sum(colSums(overlap_n > 0) > 0),
    pairs               = pairs_df
  )
}
