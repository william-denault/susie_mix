# Small history-only fixture: no model fitting, Slurm calls or production files.
code<-normalizePath("output/slider_online_015_2026-11-05",winslash="/")
source(file.path(code,"engine.R"));mb_load(code)
project<-tempfile("history_batch_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(project)
root<-file.path(project,"results_slide_prior_em");dir.create(root)
branch<-file.path(project,"workflows/slider_online_015_025");dir.create(branch,recursive=TRUE)
priors<-data.frame(tissue="Brain",matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
old<-data.frame(iteration=14L,source_iteration=13L,priors,update_method="slider_prior_active_counts_v1")
mb_save(old,file.path(branch,"original_prior_history.csv"),TRUE);mb_save(old,file.path(root,"prior_history.csv"),TRUE)
mb_save(list(source_project=project),file.path(branch,"model.rds"))
mb_save(list(stage=3L,priors=priors),file.path(branch,"state.rds"))
row<-function(iter,stage,phase,batch,role) data.frame(iteration=iter,stage=stage,phase=phase,
    batch=batch,prior_role=role,rho=if(phase=="full")1 else .3,priors)
initial<-row(15L,0L,"minibatch",0L,"input")
mb_save(initial,file.path(branch,"initial_priors.csv"),TRUE)
dirs<-c("../../results_slide_prior_em/iteration_015/batches/batch_01",
        "../../results_slide_prior_em/iteration_015/batches/batch_02",
        "../../results_slide_prior_em/iteration_021/batches/full_step")
mb_save(data.frame(stage=1:3,directory=dirs),file.path(branch,"stages.csv"),TRUE)
for(i in 1:3)dir.create(file.path(branch,dirs[i]),recursive=TRUE,showWarnings=FALSE)
for(i in c(15,21))writeLines("fixture",file.path(root,sprintf("iteration_%03d/ONLINE_OWNER.json",i)))
mb_save(initial,file.path(branch,dirs[1],"start_history.csv"),TRUE)
mb_save(list(fixture=TRUE),file.path(branch,dirs[1],"input.rds"))
mb_save(row(15L,1L,"minibatch",1L,"after_update"),file.path(branch,dirs[1],"history.csv"),TRUE)
mb_save(row(15L,2L,"minibatch",2L,"after_update"),file.path(branch,dirs[2],"history.csv"),TRUE)
mb_save(row(21L,2L,"full",0L,"input"),file.path(branch,dirs[3],"start_history.csv"),TRUE)
# Old-format full-step history is also normalized when publishing snapshots.
mb_save(row(21L,3L,"full",1L,"after_update"),file.path(branch,dirs[3],"history.csv"),TRUE)
frozen<-tools::md5sum(file.path(branch,c("state.rds","model.rds","original_prior_history.csv")))
mb_history(branch)
h<-mb_read(file.path(root,"prior_history.csv"));h14<-mb_read(file.path(root,"iteration_015/prior_history.csv"))
h20<-mb_read(file.path(root,"iteration_021/prior_history.csv"));long<-mb_read(file.path(root,"prior_history_long.csv"))
stopifnot(all(c("iteration","batch") %in% names(h14)),
 identical(h14$iteration,c(14L,15L,15L,15L)),identical(h14$batch,c(0L,0L,1L,2L)),
 all(h$batch[h$phase=="full"]==0L),all(h20$batch[h20$phase=="full"]==0L),
 all(long$batch[long$phase=="full"]==0L),
 identical(sort(h$prior_role[h$iteration==21L]),c("after_update","input")),
 all(abs(as.matrix(h[,spe_columns()])-1/17)<1e-14),
 identical(frozen,tools::md5sum(names(frozen))))
hashes<-tools::md5sum(file.path(root,c("prior_history.csv","prior_history_long.csv","iteration_015/prior_history.csv")))
mb_history(branch);stopifnot(identical(hashes,tools::md5sum(names(hashes))))
cat("PASS: iteration_015 history columns; full EM always batch 0; online batches retained;\n",
    "input/update roles distinct; root, per-scan and long histories agree; weights and fitting state unchanged; idempotent publication.\n")
