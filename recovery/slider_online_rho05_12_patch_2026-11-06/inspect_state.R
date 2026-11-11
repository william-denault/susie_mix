# Read checkpoint metadata only; do not fit, commit, publish, or modify EM files.
a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
branch <- a[1L]; destination <- a[2L]
st <- readRDS(file.path(branch, "state.rds"))
stopifnot(length(st$stage) == 1L, is.finite(st$stage), st$stage == as.integer(st$stage),
          st$stage >= 9L, is.data.frame(st$latest), nrow(st$latest) == 18263L,
          !anyNA(st$latest$gene), !anyDuplicated(st$latest$gene),
          all(nzchar(st$latest$gene)), all(file.exists(st$latest$path)),
          all(file.exists(st$latest$priors)))
write.csv(data.frame(stage = st$stage, genes = nrow(st$latest)),
          file.path(destination, "checkpoint.csv"), row.names = FALSE)
writeLines(as.character(st$latest$gene), file.path(destination, "cohort.txt"))
