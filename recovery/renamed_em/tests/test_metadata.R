args <- commandArgs(TRUE)
source(file.path(args[1],"preflight.R"))
map <- data.frame(gene=c("OLD__E1","OLD__E2"),expression_symbol="OLD",expression_id=c("E1.2","E2.4"),
  stable_id=c("E1","E2"),annotation_symbol=c("NEW1","NEW2"),chromosome="1",start=c(10,30),end=c(20,40),strand="+")
expr <- data.frame(Name=c("E1.2","E2.4"),Description="OLD")
ann <- data.frame(ensembl=c("E1","E2"),gene_name=c("NEW1","NEW2"),chromosome="chr1",start=c(10,30),end=c(20,40),strand="+")
stopifnot(all(check_metadata(map,expr,ann)$status=="PASS"))
bad <- expr;bad$Name[2]<-"E1.2"
stopifnot(all(check_metadata(map,bad,ann)$status=="FAIL"))
bad <- ann;bad$start[2]<-31
stopifnot(identical(check_metadata(map,expr,bad)$status,c("PASS","FAIL")))
bad <- expr;bad$Name[1]<-"E1.3"
stopifnot(identical(check_metadata(map,bad,ann)$status,c("FAIL","PASS")))
cat("PASS: exact stable-ID mappings, duplicate old names, missing/ambiguous rows, version and locus drift.\n")
