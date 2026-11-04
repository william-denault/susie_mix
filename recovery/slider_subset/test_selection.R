old<-new.env();new<-new.env()
source('output/slider_inplace_2026-10-29/preflight.R',local=old)
source('output/slider_inplace_646_2026-10-30/preflight.R',local=new)
original<-old$read_rename_map('output/slider_inplace_2026-10-29')
selected<-new$read_rename_map('output/slider_inplace_646_2026-10-30')
excluded<-read.csv('output/slider_inplace_646_2026-10-30/excluded_candidates.csv')
expression<-data.frame(Name=original$expression_id,Description=original$expression_symbol)
annotation<-data.frame(ensembl=original$stable_id,gene_name=original$annotation_symbol,
 chromosome=as.character(original$chromosome),start=original$start,end=original$end,strand=original$strand)
affected<-annotation$ensembl %in% excluded$stable_id
annotation$end[affected]<-annotation$end[affected]+1L
before<-old$check_metadata(original,expression,annotation)
after<-new$check_metadata(selected,expression,annotation)
stopifnot(sum(before$status=='FAIL')==34L,setequal(before$gene[before$status=='FAIL'],excluded$gene),
 nrow(after)==646L,all(after$status=='PASS'))
# Compare content after resetting row names: selected candidates are unchanged.
expected<-original[!original$gene %in% excluded$gene,];rownames(expected)<-NULL
stopifnot(identical(selected,expected))
bad<-annotation;bad$start[match(selected$stable_id[1],bad$ensembl)]<- -1L
stopifnot(new$check_metadata(selected,expression,bad)$status[1]=='FAIL')
cat('PASS: exact 646-ID subset passes; 34 endpoint mismatches are explicitly deferred; other mismatches remain blocked.\n')
