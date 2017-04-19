library("data.table")
library("GEOquery")
library("affy")
library("oligo")
library("pd.drosophila.2")
library(dplyr)


setwd("~/IMG/HP1.kd.MA.Kc/")
dir.create("CELs")
hmm.data <- read.delim("HMM.data.txt")
hmm.data2 <- fread("HMM.data.txt")
affymet <- read.delim("affymetrix_key_GPL1322-26772.txt", sep="\t", comment.char = "#")

cels <- list.files("GSE18092data/", pattern = "CEL")
sapply(paste0("GSE18092data/", cels), gunzip)
cels <- list.files("GSE18092data/", pattern = "CEL")

setwd("GSE18092data/")
raw.data.affy <- ReadAffy(cels[1:10], verbose = T)
raw.data=read.celfiles(cels[1:10])
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
