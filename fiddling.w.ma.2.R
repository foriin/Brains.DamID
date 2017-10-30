library("data.table")
library("GEOquery")
library("affy")
library("oligo")
library("pd.drosophila.2")
library(dplyr)


setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/Schwaiger_dsHP1_microarrays/ExpressionTiling/")

cels <- list.files(".", pattern = "CEL$")
wt.cel <- TileReadCel(cel.filename=cels[1],
                      bpmap.filename="Dm35b_MR_v02-3_BDGPv4h.new.bpmap", group="Dm", verbose=FALSE)
wt <- TileReadCel(cel.filename = cels[1:2],
                     bpmap.filename = "Dm35b_MR_v02-3_BDGPv4h.new.bpmap",
                     group = "", gc=F, normalize = T, verbose = F)

dsHP1 <- TileReadCel(cel.filename = cels[3:4],
                  bpmap.filename = "Dm35b_MR_v02-3_BDGPv4h.new.bpmap",
                  group = "", gc=F, normalize = T, verbose = F)
bpmap <- "Dm35b_MR_v02-3_BDGPv4h.new.bpmap"
cel <- list.files(pattern = "CEL$")
Data <- ReadAffy(widget = T)
eset <- expresso(Data, normalize.method="qspline",
                 bgcorrect.method="rma",pmcorrect.method="pmonly",
                 summary.method="liwong")
EXP <- AnalyzeTilingCelFiles(dir(pattern = ".cel|.CEL"), "Dm35b_MR_v02-3_BDGPv4h.new.bpmap")

setwd("GSE18092data/")
raw.data.affy <- ReadAffy("../Schwaiger_dsHP1_microarrays/ExpressionTiling/GSM452294_sch20070412dtr_01_ctrl_F_Kc_none_7d_97.CEL", verbose = T)
raw.data=read.celfiles("../Schwaiger_dsHP1_microarrays/ExpressionTiling/GSM452294_sch20070412dtr_01_ctrl_F_Kc_none_7d_97.CEL")
raw.data@phenoData@data

data.rma.norm=rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(data.rma.norm)

#Format values to 5 decimal places
rma=format(rma, digits=5)

lol <- affymet
affymet[18917,] <- lol[18918,]
affymet[18918,] <- lol[18917,]

rma2 <- as.data.frame(rma)

ma_data <- cbind(affymet[, c(1,2, 4, 8:11)], rma2)

ma_data[, 1:7] <- sapply(ma_data[, 1:7], as.character)

ma_data[, 8:17] <- sapply(ma_data[, 8:17], function(x) as.numeric(as.character(x)))


ma_data2 <- ma_data %>% 
  mutate(WT.exp.avg = rowMeans(cbind(.[,8:12]), na.rm = T),
         HP1KD.exp.avg = rowMeans(cbind(.[, 13:17]), na.rm = T),
         diff = WT.exp.avg - HP1KD.exp.avg) %>% 
  select(1:7, 18:20) %>% 
  arrange(diff)
