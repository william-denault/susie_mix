rm(list = ls())
source("/project2/mstephens/wdenault/susie_mix/script/sim/sim_workhorse.R")
dir.create("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/", recursive = TRUE, showWarnings = FALSE)
if (file.exists("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/additive_recessive_dominant_add2_rec1_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData")) {
  load("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/additive_recessive_dominant_add2_rec1_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData")
} else {
  results <- list()
}
start <- length(results) + 1
if (start <= 400) for (o in start:400) {
  replication <- 0 + o
  sim_seed <- 1e+06 + replication
  temp <- tryCatch(
    sim_mix(pve = 0.2, n = 500, L = 10, L_add = 2, L_rec = 1, L_dom = 1, all_additive = FALSE, seed = sim_seed, temp_dir = "/project2/mstephens/wdenault/susie_mix/temp_plink/"),
    error = function(e) list(error = conditionMessage(e))
  )
  temp$seed <- sim_seed
  temp$replication <- replication
  results[[o]] <- temp
  save(results, file = "/project2/mstephens/wdenault/susie_mix/simulation results/chunks/additive_recessive_dominant_add2_rec1_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData.tmp")
  if (!file.rename("/project2/mstephens/wdenault/susie_mix/simulation results/chunks/additive_recessive_dominant_add2_rec1_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData.tmp", "/project2/mstephens/wdenault/susie_mix/simulation results/chunks/additive_recessive_dominant_add2_rec1_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData")) stop("Could not replace the checkpoint.")
}
