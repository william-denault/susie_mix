rm(list = ls())
source("/project2/mstephens/wdenault/susie_mix/script/sim/sim_workhorse.R")
dir.create("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/", recursive = TRUE, showWarnings = FALSE)
if (file.exists("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/recessive_add0_rec2_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData")) {
  load("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/recessive_add0_rec2_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData")
} else {
  results <- list()
}
start <- length(results) + 1
if (start <= 400) for (o in start:400) {
  replication <- 0 + o
  sim_seed <- 1e+06 + replication
  temp <- tryCatch(
    sim_mix(pve = 0.1, n = 500, L = 10, L_add = 0, L_rec = 2, L_dom = 0, all_additive = FALSE, seed = sim_seed, temp_dir = "/project2/mstephens/wdenault/susie_mix/temp_plink/"),
    error = function(e) list(error = conditionMessage(e))
  )
  temp$seed <- sim_seed
  temp$replication <- replication
  results[[o]] <- temp
  save(results, file = "/project2/mstephens/wdenault/susie_mix/simulation results/chunks/recessive_add0_rec2_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData.tmp")
  if (!file.rename("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/recessive_add0_rec2_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData.tmp", "/project2/mstephens/wdenault/susie_mix/simulation results/chunks/recessive_add0_rec2_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData")) stop("Could not replace the checkpoint.")
}
