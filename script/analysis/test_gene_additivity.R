source("/project2/mstephens/wdenault/susie_mix/script/analysis/descriptive results.R", echo=TRUE)
library(data.table)

res_all <- as.data.table(res_summary)

res_all[
  , row_total_cs := n_add + n_rec + n_dom
]

gene_cs_all <- res_all[
  !is.na(row_total_cs),
  .(
    nonadditive_cs = sum(n_rec + n_dom),
    total_cs = sum(row_total_cs),
    n_tissues_tested = uniqueN(tissue),
    n_tissues_with_cs = uniqueN(tissue[row_total_cs > 0])
  ),
  by = gene
][total_cs > 0]

gene_cs_all[
  , nonadditive_ratio := nonadditive_cs / total_cs
]
gene_cs_all[
  , p_value := mapply(
    function(x, n) {
      binom.test(
        x,
        n,
        p = overall_nonadditive_ratio
      )$p.value
    },
    nonadditive_cs,
    total_cs
  )
]

gene_cs_all[
  , direction := fifelse(
    nonadditive_ratio > overall_nonadditive_ratio,
    "more non-additive",
    "more additive"
  )
]

gene_cs_all[, p_adjusted := p.adjust(p_value, method = "BH")]

gene_cs_all[order(p_adjusted)][1:100]




hist(gene_cs_all$nonadditive_ratio,nclass = 300)
