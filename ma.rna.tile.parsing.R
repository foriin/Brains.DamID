library(Starr)
library(limma)
rm(list = ls())


setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/Schwaiger_dsHP1_microarrays/ExpressionTiling/")




bpmap <- readBpmap("Dm35b_MR_v02-3_BDGPv4h.new.bpmap")
cels <- dir(pattern = "CEL$")
names <- c("Kc_control_1", "Kc_control2", "Kc_dsHP1_1", "Kc_dsHP1_2")
type <- c("CONTROL", "CONTROL", "KD", "KD")

dsHP1 <- readCelFile(bpmap, cels, names, type)
E <- dsHP1@assayData$exprs
N <- backgroundCorrect.matrix(E, method = "normexp")
dsHP1.df <- data.frame(featureData(dsHP1)@data, N)

write.table(dsHP1.df, "dsHP1.norm.exprs.csv", quote = F, row.names = F, sep = '\t', dec = ',', col.names = T)


