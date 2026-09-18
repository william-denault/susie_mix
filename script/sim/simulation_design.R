# Shared design and checkpoint naming for the five-effect simulations.
sim_schema_version <- 2L
sim_count_columns <- c("L_add", "L_rec", "L_dom", "L_prec", "L_pdom")
sim_effect_delta <- c(additive = 0, recessive = -1, dominant = 1,
                      partial_recessive = -0.5, partial_dominant = 0.5)
sim_effect_labels <- c("Additive", "Recessive", "Dominant",
                       "Partial recessive", "Partial dominant")

sim_scenario_name <- function(counts, display = FALSE) {
  stopifnot(length(counts) == 5L, all(is.finite(counts)),
            all(counts >= 0), all(counts == floor(counts)),
            sum(counts) %in% 1:5, sum(counts > 0) <= 3L)
  labels <- if (display) sim_effect_labels else names(sim_effect_delta)
  selected <- labels[counts > 0]
  if (display && length(selected) == 1L) return(paste(selected, "only"))
  if (display) selected[-1L] <- tolower(selected[-1L])
  paste(selected, collapse = if (display) " + " else "_")
}

sim_conditions <- function() {
  counts <- expand.grid(rep(list(0:5), 5), KEEP.OUT.ATTRS = FALSE)
  names(counts) <- sim_count_columns
  K <- rowSums(counts)
  counts <- counts[K %in% 1:5 & rowSums(counts > 0) <= 3L, , drop = FALSE]
  counts <- counts[do.call(order, c(list(rowSums(counts)), unname(counts))), ]
  rownames(counts) <- NULL
  data.frame(name = apply(counts, 1L, sim_scenario_name), K = rowSums(counts),
             counts, all_additive = FALSE, row.names = NULL)
}

sim_configuration <- function(counts) {
  paste0(c("add", "rec", "dom", "prec", "pdom"), counts, collapse = "_")
}

sim_checkpoint_name <- function(condition, n, L, pve, seed_base, reps, chunk) {
  paste0(condition$name, "_", sim_configuration(as.numeric(condition[sim_count_columns])),
         "_n", n, "_L", L, "_pve", pve, "_seed", seed_base,
         "_reps", reps, "_chunk", chunk, ".RData")
}

# Read the old three-count filenames too, for explicitly requested legacy plots.
sim_parse_checkpoints <- function(files) {
  pattern <- paste0("^.*_add([0-9]+)_rec([0-9]+)_dom([0-9]+)",
                    "(?:_prec([0-9]+)_pdom([0-9]+))?_n([0-9]+)_L([0-9]+)",
                    "_pve([^_]+)_seed([^_]+)_reps([0-9]+)_chunk([0-9]+)\\.RData$")
  parts <- regmatches(basename(files), regexec(pattern, basename(files), perl = TRUE))
  if (any(lengths(parts) != 12L)) stop("Unrecognized simulation checkpoint filename.")
  values <- lapply(parts, function(z) {
    z <- z[-1L]
    z[4:5][z[4:5] == ""] <- "0"
    as.numeric(z)
  })
  info <- as.data.frame(do.call(rbind, values))
  names(info) <- c(sim_count_columns, "n", "fit_L", "pve", "seed_base", "requested_reps", "chunk")
  if (anyNA(info)) stop("Invalid numeric fields in checkpoint filenames.")
  info$file <- files
  info
}
