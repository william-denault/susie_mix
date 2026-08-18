source("/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/workhorse.R", echo=TRUE)
target_gene="ZNF232"
run_susie_gene(target_gene= target_gene)
target_gene="C1QTNF9B"
grep( target_gene,lf)
