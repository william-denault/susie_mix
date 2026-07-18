# filter_protein_coding_genes.R
#
# Reads genes_all.txt (from get_all_genes.R) and restricts it to
# protein-coding genes only, using gene_biotype from the same GTF
# annotation already used elsewhere in the pipeline. This is the
# authoritative source for biotype -- filtering on the expression
# file's Description column alone can't distinguish protein-coding
# from lncRNA/pseudogene etc.
#
# Output:
#   genes_protein_coding.txt   -- final gene list to feed generate_chunks.R
#   gtf_gene_biotypes.tsv      -- full gene_id/gene_name/gene_biotype table
#   biotype_counts.txt         -- sanity-check counts per biotype

library(data.table)

datadir  <- "/project2/mstephens/gtex"
gtf_file <- file.path(datadir,
                      "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")

cat("Reading GTF (gene-level rows only take a moment)...\n")
gtf <- fread(cmd = paste("zcat",gtf_file),
             sep = "\t",header = FALSE,quote = "",
             col.names = c("seqname","source","feature","start","end",
                           "score","strand","frame","attribute"))

# GTF files can carry a few leading "#!" header lines -- drop any stray
# non-data rows just in case.
gtf <- subset(gtf,feature == "gene")
cat("Gene-level rows in GTF:",nrow(gtf),"\n")

# Attribute field looks like:
#   gene_id "ENSG00000000003"; gene_version "14"; gene_name "TSPAN6"; ... gene_biotype "protein_coding";
get_field <- function(x,field) {
  pattern <- paste0('.*',field,' "([^"]+)".*')
  ifelse(grepl(pattern,x),sub(pattern,"\\1",x),NA_character_)
}

gene_tab <- data.frame(
  gene_id      = get_field(gtf$attribute,"gene_id"),
  gene_name    = get_field(gtf$attribute,"gene_name"),
  gene_biotype = get_field(gtf$attribute,"gene_type"),
  stringsAsFactors = FALSE
)

#fwrite(gene_tab,"gtf_gene_biotypes.tsv",sep = "\t")

biotype_counts <- sort(table(gene_tab$gene_biotype),decreasing = TRUE)
cat("Top biotypes:\n")
print(head(biotype_counts,15))
#writeLines(capture.output(print(biotype_counts)),"biotype_counts.txt")

# Restrict to protein-coding.
pc <- subset(gene_tab,gene_biotype == "protein_coding")
cat("\nProtein-coding genes in GTF:",nrow(pc),"\n")
cat("Unique protein-coding symbols:",length(unique(pc$gene_name)),"\n")

# Intersect with genes_all.txt so we only keep genes that are also
# actually present in the GTEx expression file (guards against any
# annotation-version mismatch between the two files).
genes_all <- readLines("/project2/mstephens/wdenault/susie_mix/data/genes_all.txt")
genes_pc  <- sort(intersect(unique(pc$gene_name),genes_all))

cat("\nAfter intersecting with genes_all.txt:",length(genes_pc),"genes\n")
missing_from_expr <- setdiff(unique(pc$gene_name),genes_all)
if (length(missing_from_expr) > 0)
  cat(length(missing_from_expr),
      "protein-coding symbols from the GTF were NOT found in genes_all.txt",
      "(likely annotation-version mismatch -- worth spot-checking a few)\n")

writeLines(genes_pc,"/project2/mstephens/wdenault/susie_mix/data/genes_protein_coding.txt")
cat("Wrote",length(genes_pc),"genes to genes_protein_coding.txt\n")
