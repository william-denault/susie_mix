if (nzchar(Sys.getenv("COLOC_TEST_LIBRARY"))) .libPaths(c(Sys.getenv("COLOC_TEST_LIBRARY"), .libPaths()))
source("script/analysis/em_gwas_coloc_utils.R")
cg_need(c("coloc", "data.table", "jsonlite", "digest", "R.utils"))
v <- data.frame(chromosome="1",position=1:120*1000,ref="A",alt="G",build="GRCh38")
v$snp <- cg_key(v$chromosome,v$position,v$ref,v$alt)
make_qtl <- function(lead) {
  bf <- matrix(-2,1,120); bf[1,lead]<-25;colnames(bf)<-v$snp
  list(lbf=bf,alpha=exp(bf-cg_lse(bf)),variants=v,components=1L)
}
bundle <- list(schema_version=1L,metadata=list(iteration=18L),em=make_qtl(10),vanilla=make_qtl(40))
s <- list(trait_id="external_gene",target_gene="GENE",case_tissues=list("Blood"),
  study="Independent cohort",dataset_id="QTDtest",format="eqtl_catalogue_susie",trait_group="eqtl_replication",
  sample_overlap="distinct_cohort_design",bf_scale="natural_log_absolute",snp_coverage="full_finemapped_region")
g <- data.frame(component=2L,chromosome=v$chromosome,position=v$position,effect_allele=v$alt,other_allele=v$ref,
  log_bf=as.numeric(bundle$em$lbf))
gw <- cg_published_susie(g,s)
priors <- data.frame(p1=1e-4,p2=1e-4,p12=1e-5)
stopifnot(cg_changes(cg_pair(bundle,gw,priors)$pairs)$classification=="gain")
swapped<-bundle;swapped$em<-bundle$vanilla;swapped$vanilla<-bundle$em
stopifnot(cg_changes(cg_pair(swapped,gw,priors)$pairs)$classification=="loss")
# A strong external indel absent from expression must fail coverage, not disappear.
indel<-rbind(g,data.frame(component=2L,chromosome="1",position=999,effect_allele="AT",other_allele="A",log_bf=40))
bad<-cg_pair(bundle,cg_published_susie(indel,s),priors)$pairs
stopifnot(all(bad$status=="insufficient_posterior_coverage"))
err<-try(cg_published_susie(rbind(g,g[1,]),s),silent=TRUE)
stopifnot(inherits(err,"try-error"))
stopifnot(is.null(cg_case_scope(s,data.frame(gene="GENE",tissue="Blood"))),
  cg_case_scope(s,data.frame(gene="OTHER",tissue="Blood"))=="different_target_gene",
  cg_case_scope(s,data.frame(gene="GENE",tissue="Liver"))=="outside_prespecified_tissues")
gtex<-s;gtex$study<-"GTEx_v8"
stopifnot(cg_case_scope(gtex,data.frame(gene="GENE",tissue="Blood"))=="replication_independence_not_established")

job<-tempfile("molecular-job-");dir.create(job);dir.create(file.path(job,"gwas"));dir.create(file.path(job,"bundles"))
saveRDS(list(iteration=18L),file.path(job,"job_metadata.rds"))
bp<-file.path(job,"bundles/test.rds");saveRDS(bundle,bp)
cases<-data.frame(case_id=c("case1","case2"),gene=c("GENE","OTHER"),tissue="Blood",bundle="bundles/test.rds",
  md5=unname(tools::md5sum(bp)),selection_group=c("shifted","unshifted_control"))
write.csv(cases,file.path(job,"cases.csv"),row.names=FALSE)
write.csv(data.frame(case_id=cases$case_id,chromosome="1",start=1,end=200000),file.path(job,"regions.csv"),row.names=FALSE)
jsonlite::write_json(list(sources=list(s)),file.path(job,"sources.json"),auto_unbox=TRUE)
gf<-file.path(job,"gwas/external_gene.tsv.gz");data.table::fwrite(g,gf,sep="\t")
jsonlite::write_json(list(source=s,output_sha256=digest::digest(file=gf,algo="sha256")),file.path(job,"gwas/external_gene.json"),auto_unbox=TRUE)
out<-run_em_gwas_coloc(job)
audit<-read.csv(file.path(out,"audit.csv"))
change<-read.csv(file.path(out,"paired_changes.csv"))
stopifnot(nrow(audit)==2,all(audit$mode=="published"),any(audit$status=="not_in_scope"),
  any(audit$status=="completed"),all(change$classification=="gain"),all(change$selection_group=="shifted"),
  all(change$trait_group=="eqtl_replication"),all(paste0("delta_H",0:4) %in% names(change)))
cat("PASS: molecular H0-H4 gains/losses, full-region coverage, same-gene/tissue scope, GTEx exclusion and published-BF disk workflow.\n")
