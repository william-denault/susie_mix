# get_all_genes.R
#
# Extracts the full list of genes from the GTEx gene-reads file, without
# loading the (huge) expression matrix itself -- only the first two
# columns (Ensembl ID = "Name", symbol = "Description") are read.
#
# Output:
#   genes_all.txt        -- one gene symbol per line (deduplicated)
#   gene_id_map.tsv       -- Name (Ensembl ID) <-> Description (symbol), full map
#   duplicate_symbols.txt -- symbols that map to >1 Ensembl ID (needs a look)

library(data.table)

datadir <- "/project2/mstephens/gtex"

cat("Reading gene ID / symbol columns only (fast path via fread select=)...\n")
gene_info <- fread(file.path(datadir,
                             "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"),
                   sep = "\t", skip = 2, header = TRUE,
                   select = c("Name","Description"),
                   showProgress = TRUE)
class(gene_info) <- "data.frame"

cat("Total gene rows:",nrow(gene_info),"\n")
cat("Unique Ensembl IDs:",length(unique(gene_info$Name)),"\n")
cat("Unique symbols:",length(unique(gene_info$Description)),"\n")

# Write full ID <-> symbol map, for reference / debugging.
#fwrite(gene_info,"gene_id_map.tsv",sep = "\t")

# Flag symbols that map to more than one Ensembl ID -- if you loop by
# symbol (as the original script does), these are ambiguous.
dup_counts <- table(gene_info$Description)
dup_symbols <- names(dup_counts[dup_counts > 1])
if (length(dup_symbols) > 0) {
  cat(length(dup_symbols),
      "duplicate symbols found (map to >1 Ensembl ID) -- see duplicate_symbols.txt\n")
  writeLines(dup_symbols,"duplicate_symbols.txt")
} else {
  cat("No duplicate symbols found.\n")
}

# Write the deduplicated gene symbol list used to drive the analysis.
genes_all <- sort(unique(gene_info$Description))
writeLines(genes_all,"/project2/mstephens/wdenault/susie_mix/data/genes_all.txt")
cat("Wrote",length(genes_all),"unique gene symbols to genes_all.txt\n")
