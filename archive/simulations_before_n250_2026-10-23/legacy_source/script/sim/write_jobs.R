# Generate R jobs. Run this generator once; it does not run the simulations.

project_dir <- "/project2/mstephens/wdenault/susie_mix/"

workhorse <- paste0(project_dir, "script/sim/sim_workhorse.R")
temp_dir  <- paste0(project_dir, "temp_plink/")
path      <- paste0(project_dir, "script/sim/jobs/")
chunk_dir <- paste0(project_dir, "simulation results/chunks/")

pve_values <- c(0.10, 0.20, 0.30, 0.40)
n <- 500
L <- 10

chunks_per_cell <- 1
reps_per_chunk <- 400
seed_base <- 1000000

# Enumerate every allocation of 1 to 5 causal SNPs across the three codings.
# L_add, L_rec and L_dom are COUNTS. Zero excludes that generating coding.
conditions <- data.frame()

for (K in 1:5) {
  for (L_add in 0:K) {
    for (L_rec in 0:(K - L_add)) {

      L_dom <- K - L_add - L_rec

      name <- paste(
        c("additive", "recessive", "dominant")[c(L_add, L_rec, L_dom) > 0],
        collapse = "_"
      )

      conditions <- rbind(conditions, data.frame(
        name = name,
        K = K,
        L_add = L_add,
        L_rec = L_rec,
        L_dom = L_dom,
        all_additive = FALSE
      ))
    }
  }
}

print(conditions)

dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(conditions, paste0(path, "conditions.csv"), row.names = FALSE)

tt <- 1

for (i in seq_len(nrow(conditions))) {
  for (pve in pve_values) {
    for (chunk in seq_len(chunks_per_cell)) {

      output_file <- paste0(path, "sim_job_", tt, ".R")

      out_rdata <- paste0(
        chunk_dir, conditions$name[i],
        "_add", conditions$L_add[i],
        "_rec", conditions$L_rec[i],
        "_dom", conditions$L_dom[i],
        "_n", n, "_L", L, "_pve", pve,
        "_seed", seed_base, "_reps", reps_per_chunk,
        "_chunk", chunk, ".RData"
      )

      fileConn <- file(output_file, open = "w")

      writeLines("rm(list = ls())", fileConn)
      writeLines(paste0('source("', workhorse, '")'), fileConn)
      writeLines(paste0(
        'dir.create("', chunk_dir,
        '", recursive = TRUE, showWarnings = FALSE)'
      ), fileConn)

      writeLines(paste0('if (file.exists("', out_rdata, '")) {'), fileConn)
      writeLines(paste0('  load("', out_rdata, '")'), fileConn)
      writeLines("} else {", fileConn)
      writeLines("  results <- list()", fileConn)
      writeLines("}", fileConn)

      writeLines("start <- length(results) + 1", fileConn)
      writeLines(paste0(
        "if (start <= ", reps_per_chunk,
        ") for (o in start:", reps_per_chunk, ") {"
      ), fileConn)

      # Same seed across PVE values for each causal configuration.
      # Different chunks use different replication numbers.
      writeLines(paste0(
        "  replication <- ", (chunk - 1) * reps_per_chunk, " + o"
      ), fileConn)
      writeLines(paste0("  sim_seed <- ", seed_base, " + replication"), fileConn)

      # sim_mix() resets its own seed: pass seed explicitly.
      writeLines("  temp <- tryCatch(", fileConn)
      writeLines(paste0(
        "    sim_mix(pve = ", pve, ", n = ", n, ", L = ", L,
        ", L_add = ", conditions$L_add[i],
        ", L_rec = ", conditions$L_rec[i],
        ", L_dom = ", conditions$L_dom[i],
        ", all_additive = ", conditions$all_additive[i],
        ', seed = sim_seed, temp_dir = "', temp_dir, '"),'
      ), fileConn)
      writeLines("    error = function(e) list(error = conditionMessage(e))", fileConn)
      writeLines("  )", fileConn)

      # Failed simulations are saved too, without replacing them with a new draw.
      writeLines("  temp$seed <- sim_seed", fileConn)
      writeLines("  temp$replication <- replication", fileConn)
      writeLines("  results[[o]] <- temp", fileConn)
      # Finish writing the new save before replacing the last complete checkpoint.
      writeLines(paste0('  save(results, file = "', out_rdata, '.tmp")'), fileConn)
      writeLines(paste0(
        '  if (!file.rename("', out_rdata, '.tmp", "', out_rdata, '"))',
        ' stop("Could not replace the checkpoint.")'
      ), fileConn)
      writeLines("}", fileConn)

      close(fileConn)
      cat("Saved:", output_file, "\n")
      tt <- tt + 1
    }
  }
}

cat("Generated ", tt - 1, " scripts. Array range: 1-", tt - 1, "\n", sep = "")
