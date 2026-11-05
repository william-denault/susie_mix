# Extract gene-level protein-coding annotations from a GTF file using the
# same rules used to build genes_protein_coding.txt.
get_gene_annotations <- function(gene_file) {

  out <- fread(
    file = gene_file,
    sep = "\t",
    header = FALSE,
    skip = 1,
    quote = "",
    col.names = c(
      "chromosome",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attributes"
    )
  )

  # The protein-coding input list is built from all gene-level rows,
  # regardless of whether the GTF source is ensembl, havana, or
  # ensembl_havana.
  out <- out[
    feature == "gene"
  ]

  get_field <- function(x, field) {

    pattern <- paste0(
      ".*",
      field,
      ' "([^"]+)".*'
    )

    ifelse(
      grepl(pattern, x),
      sub(pattern, "\\1", x),
      NA_character_
    )
  }

  out[, ensembl := get_field(
    attributes,
    "gene_id"
  )]

  out[, gene_type := get_field(
    attributes,
    "gene_type"
  )]

  out[, gene_name := get_field(
    attributes,
    "gene_name"
  )]

  out <- out[
    gene_type == "protein_coding" &
      !is.na(gene_name) &
      nzchar(gene_name),
    .(
      chromosome,
      source,
      feature,
      start,
      end,
      strand,
      ensembl,
      gene_type,
      gene_name
    )
  ]

  as.data.frame(
    unique(out)
  )
}
