# Generate some scatterplots on raw counts data from DamID experiment

library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(GGally)

rm(list = ls())



setwd(paste0("~/IMG/Projects/",
             "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
             "BioHMM2.qn.full.PC.HMM3/"))
load("scatters.RData")

DATA <- fread("../CSV/01.Raw.Counts.csv")
DATA2 <- cbind(DATA[, 1:4],
               as.data.frame(
                 sapply(DATA[, 5:ncol(DATA)], function(col){
                   col/sum(col, na.rm = T) * 1e6
                 })
               ))

# remove some ugly outliers
DATA2$DAM.FB.m_25mkM4HT.1.all[DATA2$DAM.FB.m_25mkM4HT.1.all > 500] <- NA
DATA2$LAM.FB.m_25mkM4HT.1.all[DATA2$LAM.FB.m_25mkM4HT.1.all > 500] <- NA
DATA2$LAM.FB.m_25mkM4HT.2.all[DATA2$LAM.FB.m_25mkM4HT.2.all > 500] <- NA
DATA2$PC.FB.m_25mkM4HT.1.all[DATA2$PC.FB.m_25mkM4HT.1.all > 1000] <- NA
DATA2$PC.FB.m_25mkM4HT.2.all[DATA2$PC.FB.m_25mkM4HT.2.all > 1000] <- NA
# DATA2$HP1.FB.m_25mkM4HT.2.all[DATA2$HP1.FB.m_25mkM4HT.2.all > 600] <- NA

ScatCor <- function(data, lastcol = 4, pref){
  png.name <- paste0(pref, ".Scatter_Plots_and_Correlations.png")
  corplot <- ggpairs(
    data[, (lastcol + 1):ncol(data)],
    title = "Scatter Plots and Pearson Correlations",
    upper = list(
      continuous = wrap("cor", size = 15)),
    lower = list(
      continuous=wrap("smooth", colour="blue")
    ),
    diag = NULL) +
    theme_grey(base_size = 20)
  print(corplot)
}

ScatterPlotting <- function(dataSet) {
  for (j in unique(sub("^([^.]+\\.[^.]+)\\..*(\\d)\\.all$", "\\1", names(dataSet)[-c(1:4)], perl=T))) {
    repSet <- sort(grep(j, names(dataSet), value=T))
    if (length(repSet) > 1) {
      png(filename = paste("scatter_on_", j, ".png", sep=""), width = 1200, height = 1200)
            up.brd <- max(max(dataSet[[repSet[2]]], na.rm=T), max(dataSet[[repSet[1]]], na.rm=T))
            par(mar = c(6,5,4,2) + 0.1)
            par(cex=1.5)
            Cor.P <- round(cor(dataSet[[repSet[1]]], dataSet[[repSet[2]]], method="pearson", use="pairwise.complete.obs"), digits=2)
            Cor.S <- round(cor(dataSet[[repSet[1]]], dataSet[[repSet[2]]], method="spearman", use="pairwise.complete.obs"), digits=2)
            plot(x=dataSet[[repSet[1]]], y=dataSet[[repSet[2]]], cex.axis = 3,
                 xlim = c(0, 0.9*up.brd), ylim = c(0, 0.9*up.brd),
                 # xlab=sub("^([^.]+\\.[^.]+)\\..*(\\d)\\.all$", "\\1.\\2", repSet[1]),
                 # ylab=sub("^([^.]+\\.[^.]+)\\..*(\\d)\\.all$", "\\1.\\2", repSet[2]),
                 xlab = "",
                 ylab = "",
                 text(x= up.brd*0.15, y=up.brd*0.85, paste0("Pearson: ", Cor.P, "\nSpearman: ", Cor.S), cex = 3),
                 mgp = c(3, 2.4, 0))
            # x <- c(0, max(dataSet[[repSet[1]]], na.rm=T)); y <- c(0, max(dataSet[[repSet[2]]], na.rm=T))
            # lines(x, y, col = "red")
            # print(ggplot(dataSet, aes(dataSet[[i]], dataSet[[x]]))+geom_point(alpha=1/10, colour="red", size=4) + xlab(i) + ylab(x) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(dataSet[, i], na.rm=T), max(dataSet[, x], na.rm=T)*0.75, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw())
            rm(Cor.P)
            rm(Cor.S)
      dev.off()
    } else {
      print(paste("Skip make the Scatter Plots from", j, sep=" "))
    }
  }
}


# names(DATA) <- gsub("^([^.]+\\.[^.]+)\\..*(\\d)\\.all$", "\\1.\\2", names(DATA))
ScatterPlotting(DATA2)

