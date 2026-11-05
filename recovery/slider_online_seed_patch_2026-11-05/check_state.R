args<-commandArgs(TRUE)
branch<-normalizePath(args[1],winslash="/",mustWork=TRUE)
expected_genes<-if(length(args)>1L)as.integer(args[2])else 18263L
seed<-file.path(branch,"seed/snapshot_001")
state<-readRDS(file.path(branch,"state.rds"))
model<-readRDS(file.path(branch,"model.rds"))
metadata<-readRDS(file.path(seed,"metadata.rds"))
plan<-read.csv(file.path(branch,"stages.csv"))
first<-file.path(branch,plan$directory[1])
input<-readRDS(file.path(first,"input.rds"))
expected<-read.csv(file.path(seed,"expected.csv"),stringsAsFactors=FALSE)
stopifnot(state$stage==0L,model$source_iteration==14L,metadata$source_iteration==14L,
          nrow(plan)==59L,plan$iteration[1]==15L,plan$batch[1]==1L,plan$rho[1]==.3,
          input$stage==1L,input$iteration==15L,input$batch==1L,input$rho==.3,
          identical(input$previous_state_md5,unname(tools::md5sum(file.path(branch,"state.rds")))),
          length(model$genes)==expected_genes,identical(state$latest$gene,model$genes),
          identical(expected$gene,model$genes),all(state$latest$iteration==14L),
          identical(state$latest$md5,expected$md5),
          identical(normalizePath(state$latest$path,winslash="/",mustWork=FALSE),
                    normalizePath(file.path(seed,"results",paste0(model$genes,".rds")),winslash="/",mustWork=FALSE)),
          isTRUE(all.equal(state$S,metadata$counts/length(model$genes))),
          isTRUE(all.equal(state$priors,metadata$initial_priors)),
          isTRUE(all.equal(input$S,state$S)),isTRUE(all.equal(input$priors,state$priors)),
          identical(input$warm$gene,input$genes),
          identical(input$warm$md5,expected$md5[match(input$genes,expected$gene)]))
cat("Validated unchanged 014 seed, stage 0 state, 015 batch 1, priors and running counts.\n")
