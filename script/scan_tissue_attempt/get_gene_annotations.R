# This is a function I use to Extract the gene annotations from the GTF
# ("gene transfer format") file. Here we keep only the annotated gene
# transcripts for protein-coding genes as defined in the Ensembl/Havana
# database.
get_gene_annotations <- function (gene_file) {
  out <- fread(file = gene_file,sep = "\t",header = FALSE,skip = 1)
  class(out) <- "data.frame"
  names(out) <- c("chromosome","source","feature","start","end","score",
                  "strand","frame","attributes")
  out <- out[c("chromosome","source","feature","start","end","strand",
               "attributes")]
  out <- transform(out,
                   chromosome = factor(chromosome),
                   source     = factor(source),
                   feature    = factor(feature),
                   strand     = factor(strand))
  #out <- subset(out,
  #              source == "ensembl_havana" &
  #                feature == "transcript")

  out <- subset(
    out,
    feature == "gene"
  )
  out <-
    transform(out,
              ensembl   = sapply(strsplit(attributes,";"),
                                 function (x) substr(x[[1]],10,24)),
              gene_type = sapply(strsplit(attributes,";"),
                                 function (x) substr(x[[3]],13,nchar(x[[3]]) - 1)),
              gene_name = sapply(strsplit(attributes,";"),
                                 function (x) substr(x[[4]],13,nchar(x[[4]]) - 1)))
  out <- out[-7]
  out <- transform(out,gene_type = factor(gene_type))
  out <- subset(out,gene_type == "protein_coding")
  rownames(out) <- NULL
  return(out)
}
