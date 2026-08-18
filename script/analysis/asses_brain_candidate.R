ad_gwas_tier1_genes <- unique(c(# from consensus GWAS
  "AGRN",
  "ADAM17",
  "EIF4G3",
  "PRKD3",
  "IDUA",
  "DGKQ",
  "ANKH",
  "OTULIN",

  # HLA locus
  "HLA-DQA1",
  "HLA-DRB1",
  "HLA-DRB5",

  "CLNK",
  "RHOH",
  "ADGRL3",
  "TREM2",
  "CD2AP",
  "SORT1",
  "COX7C",
  "NCK2",
  "BIN1",
  "MGAT5",
  "HBEGF",
  "HS3ST5",
  "PPP2R3A",
  "ADAMTS4",
  "MME",

  # WDR12/ICA1L/CYP20A1 locus
  "WDR12",
  "ICA1L",
  "CYP20A1",

  "TNIP1",
  "FAM193B",
  "RASGEF1C",
  "PTPRC",
  "CR1",
  "INPP5D",
  "CTSB",

  # ECHDC3/USP6NL locus
  "ECHDC3",
  "USP6NL",

  # CLU/PTK2B locus
  "CLU",
  "PTK2B",

  "TMEM184A",

  # ICA1/UMAD1 locus
  "ICA1",
  "UMAD1",

  "TMEM106B",
  "JAZF1",
  "NME8",

  # CELF1/SPI1 locus
  "CELF1",
  "SPI1",

  # MS4A locus
  "MS4A4A",
  "MS4A6A",

  "IPMK",
  "ANK3",
  "TSPAN14",
  "PICALM",
  "ABCA1",

  # ZCWPW1/NYAP1 locus
  "ZCWPW1",
  "NYAP1",

  "BLNK",
  "PLEKHA1",
  "SORL1",

  # TPCN1/RITA1/IQCD locus
  "TPCN1",
  "RITA1",
  "IQCD",

  "DOCK4",
  "EPHA1",
  "SHARPIN",

  # WDR81/SERPINF2 locus
  "WDR81",
  "SERPINF2",

  # SCIMP/RABEP1 locus
  "SCIMP",
  "RABEP1",

  # MYO15A/TOM1L2 locus
  "MYO15A",
  "TOM1L2",

  "IQCK",
  "UBFD1",
  "DOC2A",
  "FERMT2",

  # IL34/MTSS2 locus
  "IL34",
  "MTSS2",

  "GRN",
  "MAPT",
  "ABI3",
  "TSPOAP1",
  "ACE",
  "ATP8B4",
  "ADAM10",
  "APH1B",

  # SNX1/CIAO2A locus
  "SNX1",
  "CIAO2A",

  "CTSH",

  # SLC24A4/RIN3 locus
  "SLC24A4",
  "RIN3",

  "MAF",
  "PLCG2",
  "PRDM7",
  "RBCK1",
  "ABCA7",
  "VMAC",

  # KLF16/REXO1 locus
  "KLF16",
  "REXO1",

  "VAV1",
  "LRRC25",
  "SRC",
  "APP",
  "ADAMTS1",
  "CEP89",
  "CASS4",
  "APOE",
  "SIGLEC11",
  "CD33",
  "LILRB2",

  # LILRB1/LILRB4 locus
  "LILRB1",
  "LILRB4"
))


o=1
gene_name= ad_gwas_tier1_genes [o]
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)

par(mfrow=c(3,2))

k=4
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=5
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))

abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=14
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))

abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


print(gene_name)
o=o+1
#brain
#APOE
"PILRA"
"MGAT5"
"CR1"
"PRDM7"
"SIGLEC11"#in never
