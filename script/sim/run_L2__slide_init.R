project_dir <- "/project2/mstephens/wdenault/susie_mix"
Sys.setenv(SUSIE_MIX_PROJECT_DIR = project_dir)
source(file.path(project_dir, "script/sim/sim_additive_slide_init.R"))

results <- run_additive_initialization_experiment()
