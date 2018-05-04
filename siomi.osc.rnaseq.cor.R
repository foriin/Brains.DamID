library(dplyr)
library(data.table)

siomi.osc <- fread("~/IMG/Projects/HP1.ovaries/rnaseq/siomi_OSC_RNA-seq/GSM1142844_EGFP_KD_transcript.txt")

siomi.x.seq <- merge(siomi.osc, all.tpm[,1:8], by.x = "gene_id", by.y = "gene_name")

cor(siomi.x.seq$TPM.av.WT, siomi.x.seq$`EGFP KD RPKM`, use = "complete.obs")
sum(complete.cases(siomi.x.seq))
