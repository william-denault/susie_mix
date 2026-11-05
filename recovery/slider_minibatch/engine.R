# Online EM orchestration only. Scientific fitting/count/validation functions
# below are vendored unchanged from the previously tested slider bundle.
mb_load <- function(code) {
  for (f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R",
              "slide_prior_parallel.R","workhorse_slide_prior_em.R"))
    sys.source(file.path(code,"script/scan_tissue_attempt",f),globalenv())
}
mb_md5 <- function(f) unname(tools::md5sum(f))
mb_read <- function(f) read.csv(f, stringsAsFactors=FALSE)
mb_save <- function(x,f,csv=FALSE) em_atomic_write(x,f,csv=csv)
mb_weights <- function(S, previous) {
  stopifnot(is.matrix(S),identical(dimnames(S),list(previous$tissue,spe_columns())))
  w <- t(vapply(seq_len(nrow(previous)),function(i)
    spe_update_weights(S[i,],spe_prior(previous,previous$tissue[i])),numeric(17)))
  colnames(w)<-spe_columns()
  data.frame(tissue=previous$tissue,w,check.names=FALSE)
}
mb_update <- function(S, counts, batch_n, previous, rho) {
  if(length(batch_n)!=1L || batch_n<1 || !is.finite(rho) || rho<=0 || rho>1 ||
     !identical(dimnames(S),dimnames(counts)) || any(!is.finite(counts)) || any(counts<0))
    stop("Invalid online sufficient-statistic update")
  nextS <- (1-rho)*S + rho*counts/batch_n
  list(S=nextS,priors=mb_weights(nextS,previous))
}
mb_packages <- function() {
  pkgs<-c("susieR","susieRSlidePrior","data.table","matrixStats")
  setNames(vapply(pkgs,function(p) as.character(utils::packageVersion(p)),character(1)),pkgs)
}
mb_capture <- function(code, source_project, dest, expected_genes=18263L, expected_manifest=19114L) {
  # Private source runtime remains read-only; its own validator certifies the
  # last completed checkpoint, not merely the largest directory/file count.
  rt<-file.path(source_project,"slider_inplace_runtime_2026-10-29")
  e<-new.env(parent=globalenv())
  for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R",
             "slide_prior_parallel.R","cohort.R","workhorse_slide_prior_em.R"))
    sys.source(file.path(rt,"script/scan_tissue_attempt",f),e)
  normalize_body<-function(f) gsub("slider_inplace_runtime_2026-10-29/","",
                                 paste(deparse(body(f)),collapse="\n"),fixed=TRUE)
  for(f in c("spe_settings","spe_fit_counts","spe_update_weights","spe_summarize_gene","sre_validate"))
    if(!identical(normalize_body(get(f,e)),normalize_body(get(f,globalenv()))))
      stop("Source scientific helper differs from this bundle: ",f)
  if(!identical(normalize_body(e$.inplace_run_gene),normalize_body(run_slide_prior_gene)))
    stop("Source fitting function differs from this bundle")
  dirs<-list.dirs(file.path(source_project,"results_slide_prior_em"),recursive=FALSE,full.names=TRUE)
  dirs<-dirs[grepl("^iteration_[0-9]+$",basename(dirs))]
  dirs<-dirs[order(as.integer(sub("iteration_","",basename(dirs))),decreasing=TRUE)]
  dirs<-dirs[file.exists(file.path(dirs,"recovery/parallel/COMPLETE.rds")) &
             !file.exists(file.path(dirs,"storage_retired.json"))]
  if(!length(dirs)) stop("No retained, certified slider checkpoint")
  d<-dirs[1]; ctx<-e$sre_context(source_project,d,check_package=TRUE)
  if(ctx$iteration<7L) stop("Need a completed expanded-cohort checkpoint (007 or later)")
  pooled<-tryCatch(e$spp_require_complete(source_project,d),error=function(err) {
    if(grepl("Parallel completion is missing or changed",conditionMessage(err),fixed=TRUE))
      stop(structure(list(message="Source completion certificate is being refreshed; retry capture later.",call=NULL),
                     class=c("mb_snapshot_busy","error","condition")))
    stop(err)
  })
  cert<-readRDS(file.path(d,"recovery/parallel/COMPLETE.rds"))
  index<-readRDS(file.path(d,"recovery/parallel/index.rds"))
  genes<-sort(setdiff(ctx$manifest$gene,ctx$scope$excluded))
  if(length(genes)!=expected_genes || nrow(ctx$manifest)!=expected_manifest) stop("Unexpected expanded gene universe")
  provenance<-e$si_info(source_project)
  map<-mb_read(file.path(rt,"renamed_gene_map.csv"))
  if(nrow(map)!=646L || !all(map$gene %in% genes) ||
     mb_md5(file.path(rt,"renamed_gene_map.csv"))!=mb_md5(file.path(code,"renamed_gene_map.csv")))
    stop("Renamed-gene selection differs from the frozen 646-ID map")
  summaries<-index$entries[genes]
  expectations<-lapply(seq_along(genes),function(j) {
    g<-genes[j]; x<-summaries[[j]]; k<-match(g,ctx$manifest$gene)
    stamp<-cert$result_stamps[k,,drop=FALSE];rownames(stamp)<-NULL
    if(is.null(x) || !e$sre_ready(x$row$status) || !identical(x$key,cert$key) ||
       !identical(x$result_stamp,stamp) || !grepl("^[a-f0-9]{32}$",x$row$result_md5))
      stop("Uncertified source summary: ",g)
    data.frame(gene=g,source=file.path(d,"results",paste0(g,".rds")),
               size=stamp$size,mtime=stamp$mtime,md5=x$row$result_md5)
  })
  counts<-as.matrix(pooled$counts[,spe_columns()]);rownames(counts)<-pooled$counts$tissue
  counts<-counts[ctx$priors$tissue,,drop=FALSE]
  # Verify that the cached summaries and certified pooled counts agree.
  total<-Reduce(`+`,lapply(summaries,function(x) x$summary$counts))
  if(max(abs(total-counts))>1e-7) stop("Cached source counts disagree with certified pooled counts")
  metadata<-list(schema=1L,source_project=source_project,source_dir=d,
    source_iteration=ctx$iteration,settings=ctx$settings,packages=mb_packages(),
    genes=genes,manifest=ctx$manifest,used_priors=ctx$priors,
    counts=counts,initial_priors=mb_weights(counts,ctx$priors),datadir=provenance$datadir,
    cohorts=setNames(lapply(summaries,`[[`,"cohort"),genes),certificate=cert,
    initial_fits=sum(pooled$audit$issue=="fit"),map=map)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  mb_save(metadata,file.path(dest,"metadata.rds"))
  mb_save(ctx$priors,file.path(dest,"source_priors.csv"),TRUE)
  mb_save(do.call(rbind,expectations),file.path(dest,"expected.csv"),TRUE)
  mb_save(data.frame(gene=genes),file.path(dest,"genes.csv"),TRUE)
  # Copy the self-contained prior-history snapshot, not the moving root file.
  if(!file.copy(file.path(d,"prior_history.csv"),file.path(dest,"source_prior_history.csv"),overwrite=TRUE))
    stop("Could not preserve the source prior history")
  # Last write is the capture transaction's completion marker.
  mb_save(data.frame(iteration=ctx$iteration,genes=length(genes),tissue_fits=metadata$initial_fits),
          file.path(dest,"capture.csv"),TRUE)
  cat("Selected validated source",ctx$iteration,"with",length(genes),"genes and",metadata$initial_fits,"tissue fits.\n")
}
mb_initialize <- function(branch,seed) {
  statefile<-file.path(branch,"state.rds")
  if(file.exists(statefile)) return(invisible(NULL))
  m<-readRDS(file.path(seed,"metadata.rds")); expected<-mb_read(file.path(seed,"expected.csv"))
  if(!identical(m$packages,mb_packages())) stop("R package versions changed")
  stopifnot(identical(expected$gene,m$genes))
  paths<-file.path(seed,"results",paste0(m$genes,".rds"));st<-file.info(paths)
  if(anyNA(st$size) || any(st$size!=expected$size) || any(abs(as.numeric(st$mtime)-expected$mtime)>0.001))
    stop("Pinned source files are incomplete or changed")
  cfg<-m[c("settings","packages","genes","datadir","cohorts","map","source_iteration","initial_fits")]
  mb_save(cfg,file.path(branch,"model.rds"))
  S<-m$counts/length(m$genes)
  latest<-data.frame(gene=m$genes,path=paths,md5=expected$md5,
    priors=file.path(seed,"source_priors.csv"),iteration=m$source_iteration,
    stringsAsFactors=FALSE)
  st<-list(stage=0L,S=S,priors=m$initial_priors,latest=latest,seed=seed)
  mb_save(data.frame(stage=0L,phase="seed",pass=0L,batch=0L,rho=NA_real_,n_genes=length(m$genes),
    source_iteration=m$source_iteration,m$initial_priors),file.path(branch,"initial_priors.csv"),TRUE)
  mb_save(st,statefile)
  cat("Initialized running counts as certified full-data counts /",length(m$genes),".\n")
}
mb_prepare <- function(branch) {
  st<-readRDS(file.path(branch,"state.rds"));plan<-mb_read(file.path(branch,"stages.csv"))
  if(st$stage>=nrow(plan)) return(invisible(NULL))
  row<-plan[plan$stage==st$stage+1L,,drop=FALSE]
  d<-file.path(branch,row$directory); f<-file.path(d,"input.rds")
  if(file.exists(f)) {
    x<-readRDS(f)
    if(!identical(x$previous_state_md5,mb_md5(file.path(branch,"state.rds")))) stop("Stage input belongs to another state")
    return(invisible(NULL))
  }
  cfg<-readRDS(file.path(branch,"model.rds"))
  genes<-readLines(file.path(branch,row$genes_file),warn=FALSE)
  if(!length(genes)||anyDuplicated(genes)||any(!genes %in% cfg$genes)) stop("Invalid stage gene assignments")
  w<-st$latest[match(genes,st$latest$gene),,drop=FALSE]
  if(any(!file.exists(w$path))) stop("Required warm starts missing")
  x<-list(stage=as.integer(row$stage),phase=row$phase,pass=row$pass,batch=row$batch,
    rho=row$rho,genes=genes,warm=w,priors=st$priors,S=st$S,
    cohorts=cfg$cohorts[genes],previous_state_md5=mb_md5(file.path(branch,"state.rds")))
  for(sub in c("results","summaries","success","attempts","logs")) dir.create(file.path(d,sub),recursive=TRUE,showWarnings=FALSE)
  mb_save(st$priors,file.path(d,"priors.csv"),TRUE)
  mb_save(data.frame(gene=genes),file.path(d,"manifest.csv"),TRUE)
  mb_save(x,f)
  cat("Prepared stage",x$stage,x$phase,"pass",x$pass,"batch",x$batch,"genes",length(genes),"rho",x$rho,"\n")
}
mb_worker <- function(branch,d,gene,run_gene=run_slide_prior_gene) {
  x<-readRDS(file.path(d,"input.rds"));cfg<-readRDS(file.path(branch,"model.rds"))
  if(!gene %in% x$genes || !identical(cfg$packages,mb_packages())) stop("Wrong gene or changed R packages")
  settings<-cfg$settings
  if(!identical(spe_settings(),settings[names(spe_settings())])) stop("Frozen model settings differ")
  ctx<-list(iteration=x$stage,priors=x$priors,settings=settings)
  outfile<-file.path(d,"results",paste0(gene,".rds"))
  if(!file.exists(outfile)) {
    w<-x$warm[x$warm$gene==gene,,drop=FALSE]
    if(nrow(w)!=1L || mb_md5(w$path)!=w$md5) stop("Warm-start checksum changed: ",gene)
    old<-readRDS(w$path)
    sre_validate(old,list(iteration=w$iteration,priors=mb_read(w$priors),settings=settings),x$cohorts[[gene]])
    settings$max_iter<-max(5000L,settings$max_iter)
    if(!is.null(old$error) || !length(old)) out<-old else {
      map<-cfg$map[cfg$map$gene==gene,,drop=FALSE]
      args<-c(list(target_gene=if(nrow(map)) map$expression_symbol else gene,
        tissue_priors=x$priors,project_dir=file.path(branch,"code"),datadir=cfg$datadir,
        previous_result=old,temp_dir=file.path(d,"temp_plink",gene)),
        settings[setdiff(names(settings),c("schema","package_version","cohort_source"))])
      if(nrow(map)) args$target_gene_id<-map$stable_id
      out<-do.call(run_gene,args)
      if(nrow(map)) attr(out,"recovery_identity")<-as.list(map)
    }
    attr(out,"em_iteration")<-x$stage;attr(out,"fit_max_iter")<-settings$max_iter
    sre_validate(out,ctx,x$cohorts[[gene]],require_converged=TRUE)
    if(is.null(out$error)) for(t in names(out)) out[[t]]$fit_slide_prior$call<-NULL
    mb_save(out,outfile)
    rm(out,old);if(exists("args",inherits=FALSE))rm(args);gc(verbose=FALSE)
  }
  # A crash after saving a valid output is recovered here without refitting.
  out<-readRDS(outfile)
  sre_validate(out,ctx,x$cohorts[[gene]],require_converged=TRUE)
  sm<-spe_summarize_gene(out,gene,x$priors)
  sf<-file.path(d,"summaries",paste0(gene,".rds"))
  mb_save(list(gene=gene,stage=x$stage,input_md5=mb_md5(file.path(d,"input.rds")),
    counts=sm$counts,audit=sm$audit,cohort=spp_cohort(out)),sf)
  mb_save(data.frame(gene=gene,stage=x$stage,input_md5=mb_md5(file.path(d,"input.rds")),
    result_md5=mb_md5(outfile),summary_md5=mb_md5(sf),n_fits=sum(sm$audit$issue=="fit")),
    file.path(d,"success",paste0(gene,".csv")),TRUE)
}
mb_commit <- function(branch,d) {
  x<-readRDS(file.path(d,"input.rds"));st<-readRDS(file.path(branch,"state.rds"))
  # Commit payload precedes the state pointer. Repeating either operation is safe.
  f<-file.path(d,"commit.rds")
  if(file.exists(f)) {
    commit<-readRDS(f)
    if(st$stage==x$stage-1L) {
      if(!identical(mb_md5(file.path(branch,"state.rds")),x$previous_state_md5)) stop("State changed before commit recovery")
      mb_save(commit$state,file.path(branch,"state.rds"))
    } else if(st$stage!=x$stage) stop("Unexpected committed-stage state")
    return(invisible(commit))
  }
  if(st$stage!=x$stage-1L || !identical(mb_md5(file.path(branch,"state.rds")),x$previous_state_md5))
    stop("Invalid commit predecessor")
  cfg<-readRDS(file.path(branch,"model.rds"))
  genes<-sort(x$genes)
  input_md5<-mb_md5(file.path(d,"input.rds"))
  summaries<-lapply(genes,function(g) {
    rec<-mb_read(file.path(d,"success",paste0(g,".csv")))
    sf<-file.path(d,"summaries",paste0(g,".rds"));sm<-readRDS(sf)
    if(nrow(rec)!=1L || rec$gene!=g || rec$stage!=x$stage || rec$input_md5!=input_md5 ||
       rec$summary_md5!=mb_md5(sf) || sm$input_md5!=input_md5 || sm$gene!=g || sm$stage!=x$stage ||
       !identical(sm$cohort,x$cohorts[[g]])) stop("Missing/changed validated summary: ",g)
    sm
  })
  counts<-Reduce(`+`,lapply(summaries,`[[`,"counts"))
  # Full steps use exactly the established full-data M-step, including its Q check.
  if(x$phase=="full") {
    if(!setequal(genes,cfg$genes)) stop("Full EM step has an incomplete cohort")
    pooled<-spe_pool_summaries(summaries,x$priors)
    updated<-list(S=counts/length(cfg$genes),priors=pooled$priors[,c("tissue",spe_columns())])
  } else updated<-mb_update(x$S,counts,length(genes),x$priors,x$rho)
  audit<-do.call(rbind,lapply(summaries,`[[`,"audit"))
  for(g in genes) {
    j<-match(g,st$latest$gene);rec<-mb_read(file.path(d,"success",paste0(g,".csv")))
    st$latest[j,c("path","md5","priors","iteration")]<-list(
      file.path(d,"results",paste0(g,".rds")),rec$result_md5,file.path(d,"priors.csv"),x$stage)
  }
  st$stage<-x$stage;st$S<-updated$S;st$priors<-updated$priors
  hist<-data.frame(stage=x$stage,phase=x$phase,pass=x$pass,batch=x$batch,rho=x$rho,
    n_genes=length(genes),source_iteration=cfg$source_iteration,updated$priors,
    prior_max_change=apply(abs(as.matrix(updated$priors[,spe_columns()])-as.matrix(x$priors[,spe_columns()])),1,max))
  mb_save(data.frame(tissue=rownames(counts),counts),file.path(d,"counts.csv"),TRUE)
  mb_save(updated$S,file.path(d,"running_counts_per_gene.rds"))
  mb_save(updated$priors,file.path(d,"prior_after.csv"),TRUE)
  mb_save(audit,file.path(d,"audit.csv"),TRUE)
  mb_save(hist,file.path(d,"history.csv"),TRUE)
  commit<-list(stage=x$stage,state=st,history=hist,genes=genes,
    tissue_fits=sum(audit$issue=="fit"),input_md5=input_md5,
    completed_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  mb_save(commit,f)
  mb_save(st,file.path(branch,"state.rds"))
  cat("Committed stage",x$stage,"genes",length(genes),"tissue fits",commit$tissue_fits,"\n")
  invisible(commit)
}
mb_history <- function(branch) {
  st<-readRDS(file.path(branch,"state.rds"));plan<-mb_read(file.path(branch,"stages.csv"))
  initial<-mb_read(file.path(branch,"initial_priors.csv"));initial$prior_max_change<-0
  rows<-lapply(plan$directory[plan$stage<=st$stage],function(d) mb_read(file.path(branch,d,"history.csv")))
  mb_save(do.call(rbind,c(list(initial),rows)),file.path(branch,"prior_history.csv"),TRUE)
  mb_save(st$priors,file.path(branch,"latest_priors.csv"),TRUE)
  mb_save(data.frame(stage=st$stage),file.path(branch,"progress.csv"),TRUE)
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE);code<-normalizePath(a[2],winslash="/",mustWork=TRUE);mb_load(code)
  if(a[1]=="capture") tryCatch(mb_capture(code,normalizePath(a[3],winslash="/"),a[4]),
    mb_snapshot_busy=function(e) { message(conditionMessage(e));quit(status=75L) })
  else if(a[1]=="initialize") mb_initialize(a[3],a[4])
  else if(a[1]=="prepare") mb_prepare(a[3])
  else if(a[1]=="worker") mb_worker(a[3],a[4],a[5])
  else if(a[1]=="commit") mb_commit(a[3],a[4])
  else if(a[1]=="history") mb_history(a[3])
  else stop("Unknown minibatch engine mode")
}
