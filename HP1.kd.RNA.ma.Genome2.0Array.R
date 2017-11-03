# Analysis of Schwaiger HP1 KD RNA microarray data (GSE18092)


library("data.table")
library("affy")
library("oligo")
library("pd.drosophila.2")
library(tidyr)
library(dplyr)
library(limma)

rm(list=ls())
setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/HP1.kd.MA.Kc/")
load("HP1.kd.ma.RData")

# dir.create("CELs")
hmm.data <- read.delim("HMM.data.txt")
hmm.data2 <- fread("HMM.data.txt")
affymet <- read.delim("affymetrix_key_GPL1322-26772.txt", sep="\t", comment.char = "#")
euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

cels <- list.files("GSE18092data/", pattern = "CEL")
sapply(paste0("GSE18092data/", cels), gunzip)
cels <- list.files("GSE18092data/", pattern = "CEL")

setwd("GSE18092data/")
# raw.data=read.celfiles(cels[1:10])
# raw.data@phenoData@data
# data.rma.norm=rma(raw.data)

# Using these parameters I got the same results as Schwaiger 

# Raw CELs as affybatch
ma.batch <- read.affybatch(paste0("GSE18092data/", cels[1:10]))
# GC-RMA normalisation
gcrma.ma <- gcrma(ma.batch)
# MAS5 absolute detection
mas5.ma <- mas5calls(ma.batch)


mas5.ma.pres <- as.data.frame(mas5.ma@assayData$exprs) %>% transmute_all(funs(ifelse(. == "P", 1, 0)))
mas5.ma.pres.tot <- data.frame(wt.p = rowSums(mas5.ma.pres[, 1:5]), hp1.p = rowSums(mas5.ma.pres[, 6:10]))
ma.cutoff <- which(mas5.ma.pres.tot$wt.p > 2 & mas5.ma.pres.tot$hp1.p > 2)
gcrma.ma.filt <- gcrma.ma[ma.cutoff, ]

# Finding differentialy expressed genes with limma

pData(gcrma.ma.filt)
# compose design matrix
conds <- c(rep("wt", 5), rep("dsHP1", 5))
design <- model.matrix(~factor(conds, levels = c("wt", "dsHP1")))
colnames(design) <- c("wt", "dsHP1.vs.wt")
# apply regression models
fit <- lmFit(gcrma.ma, design)
fit <- eBayes(fit)
options(scipen = 9)
# tail(topTable(fit, coef = 1, n = 5000, adjust = "BH", genelist = affymet$Gene.Symbol))

# Get genes, differentially expressed in dsHP1
rma.diff <- topTable(fit, n = 20000, coef=2, genelist = affymet$ID, p.value = 0.05)

# ANALYSIS OF GENES EXPRESSION IN RELATION TO DAMID DOMAINS

#Get the important stuff out of the RMA data - the expression estimates for each array
expr = exprs(gcrma.ma.filt)

#Format values to 5 decimal places
expr = format(expr, digits=5) %>% as.data.frame() %>% mutate(ID = rownames(expr))


ma_data <- merge(affymet[, c(1,2, 4, 8:11)], expr, by = "ID")
ma_data[, 1:7] <- sapply(ma_data[, 1:7], as.character)
ma_data[, 8:17] <- sapply(ma_data[, 8:17], function(x) as.numeric(as.character(x)))


ma_data2 <- ma_data %>% 
  mutate(WT.exp.avg = log2(rowMeans(cbind(2^.[,8:12]), na.rm = T)),
         HP1KD.exp.avg = log2(rowMeans(cbind(2^.[, 13:17]), na.rm = T)),
         diff = WT.exp.avg - HP1KD.exp.avg) %>% 
  dplyr::select(1:7, 18:20) %>% 
  arrange(diff)

gene.coor <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf")

# Add gene coordinates to ma data
gene.coor <- gene.coor %>% 
  dplyr::select(1,2,4,5,7,9) %>% 
  setNames(c("chr", "source", "start", "end", "strand", "mess")) %>% 
  separate(mess, c("FBgn", "trash", "gene_name"), ";") %>% 
  mutate(FBgn = sub('.*"(.*)"', "\\1", FBgn), gene_name = sub('.*"(.*)"', "\\1", gene_name)) %>% 
  dplyr::select(-trash)

chip.key.r5 <- merge(affymet, gene.coor, by.x = "CLONE_ID_LIST", by.y = "FBgn") %>% dplyr::select(2,17, 19,20,21)

gse18092.r5 <- merge(chip.key.r5, ma_data2, by = "ID")  %>% 
  mutate(TSS = ifelse(strand == "+", start, end),
         truediff = ifelse(ID %in% rma.diff$ID, 1, 0),
         twofold = ifelse(2^WT.exp.avg / 2^HP1KD.exp.avg >= 2 | 2^WT.exp.avg / 2^HP1KD.exp.avg <= 0.5, 1, 0))

# check which genes are in domains

kc.damid.hmm <- fread("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/Kc167_multi_DamID_Filion_2010/Kc167_HP1.LAM.PC.DamID.HMM.csv")
prot.gr <- lapply(c("HP1", "LAM", "PC"), function(pr){
  GRanges(seqnames = Rle(kc.damid.hmm$chromosome[kc.damid.hmm[[pr]] == 1]),
                  ranges = IRanges(start = kc.damid.hmm$start[kc.damid.hmm[[pr]] == 1],
                                   end = kc.damid.hmm$end[kc.damid.hmm[[pr]] == 1]))
})

names(prot.gr) <- c("HP1", "LAM", "PC")


# GRanges object with TSS coordinates
tss.gr <- GRanges(seqnames = Rle(paste0("chr", gse18092.r5$chr)),
                  ranges = IRanges(start = gse18092.r5$TSS,
                                   width = 1,
                                   names = gse18092.r5$ID))

# GRanges object with all genes coordinates
genes.gr <- GRanges(seqnames = Rle(paste0("chr", gse18092.r5$chr)),
                    ranges = IRanges(start = gse18092.r5$start,
                                     end = gse18092.r5$end,
                                     names = gse18092.r5$ID))


gse18092.r5.damid.tss <- cbind(gse18092.r5, sapply(names(prot.gr), function(grname){
  gr <- prot.gr[[grname]]
  in_domain <- subsetByOverlaps(tss.gr, gr)@ranges@NAMES
  ifelse(gse18092.r5$ID %in% in_domain, 1, 0)
})) %>% mutate(chr = paste0("chr", chr), nondom = ifelse(HP1 + LAM + PC == 0, 1, 0)) %>% filter(chr %in% euc)



gse18092.r5.damid.wh.genes <- cbind(gse18092.r5, sapply(names(prot.gr), function(grname){
  gr <- prot.gr[[grname]]
  in_domain <- subsetByOverlaps(genes.gr, gr)@ranges@NAMES
  ifelse(gse18092.r5$ID %in% in_domain, 1, 0)
})) %>% mutate(chr = paste0("chr", chr), nondom = ifelse(HP1 + LAM + PC == 0, 1, 0)) %>% filter(chr %in% euc)

rma.diff.wh.genes.damid <- gse18092.r5.damid.wh.genes %>% filter(truediff == 1) %>% dplyr::select(-truediff)  %>% filter(chr %in% euc)

write.table(rma.diff.wh.genes.damid, "diff.genes.whole.domains.euc.csv", quote = F, dec = ",",
            sep = ";", row.names = F)

rma.diff.wh.genes.damid  %>% filter(chr %in% euc) %>% group_by(chr) %>% summarise(n = n())


# FUN BEGINS!


# 1. TSS

# 1.1. All genes

pdf("all.TSS.in.domains.wt.hp1kd.pdf", width = 12)
op <- par(mfrow = c(1, 4))
  # HP1
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in HP1 domains", ouline = F)
  text(1.5, 12, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1]),
                      "genes\np =",
                      format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 1],
                                         gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1],
                                         alt = "l")$p.value, digits = 2)))
  # Lamin
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$LAM == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in LAM domains", outline = F)
  text(1.5, 11.2, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1]),
                        "genes\np =",
                        format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$LAM == 1],
                                           gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1],
                                           alt = "l")$p.value, digits = 2)))
  # Polycomb
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$PC == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in PC domains", outline = F)
  text(1.5, 11, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1]),
                       "genes\np =",
                       format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$PC == 1],
                                          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1],
                                          alt = "l")$p.value, digits = 2)))
  # Conservative non-domains (WRONG, see diff genes)
  # boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$nondom == 1],
  #         gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$nondom == 1],
  #         names = c("WT", "HP1 KD"),
  #         main = "Gene expression in\nconservative non-domains", outline = F)
  # text(1.5, 5.5, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$nondom == 1]),
  #                      "genes\np =",
  #                      format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$nondom == 1],
  #                                         gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$nondom == 1],
  #                                         alt = "g")$p.value, digits = 2)))
  par(op)
dev.off()


# 1.2. Differentially expressed genes

# unlog

gse18092.r5.damid.tss <- gse18092.r5.damid.tss %>% mutate(WT.exp.avg = 2^WT.exp.avg,
                                                          HP1KD.exp.avg = 2^HP1KD.exp.avg) %>%
  filter(!(ID %in% trash_id))


pdf("dif.TSS.in.domains.wt.hp1kd.unlog.pdf", width = 12)
  op <- par(mfrow = c(1, 4))
  # HP1
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Differentially expressing\ngenes in HP1 domains", outline = F)
  text(1.5, 1800, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1]),
                        "genes\np =",
                        format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1],
                                           gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1],
                                           alt = "l")$p.value, digits = 2)))
  # LAM
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$LAM == 1 & gse18092.r5.damid.tss$truediff == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1 & gse18092.r5.damid.tss$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Differentially expressing\ngenes in LAM domains", outline = F)
  text(1.5, 1200, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1 & gse18092.r5.damid.tss$truediff == 1]),
                        "genes\np =",
                        format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$LAM == 1 & gse18092.r5.damid.tss$truediff == 1],
                                    gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$LAM == 1 & gse18092.r5.damid.tss$truediff == 1],
                                    alt = "l")$p.value, digits = 2)))
  # Polycomb
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$PC == 1 & gse18092.r5.damid.tss$truediff == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1 & gse18092.r5.damid.tss$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Differentially expressing\ngenes in PC domains", outline = F)
  text(1.5, 1000, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1 & gse18092.r5.damid.tss$truediff == 1]),
                        "genes\np =",
                        format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$PC == 1 & gse18092.r5.damid.tss$truediff == 1],
                                           gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$PC == 1 & gse18092.r5.damid.tss$truediff == 1],
                                           alt = "l")$p.value, digits = 2)))
  # Conservative non-domains
  boxplot(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 0 &
                                             gse18092.r5.damid.tss$LAM == 0 &
                                             gse18092.r5.damid.tss$PC == 0 & 
                                             gse18092.r5.damid.tss$truediff == 1],
          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 0 &
                                                gse18092.r5.damid.tss$LAM == 0 &
                                                gse18092.r5.damid.tss$PC == 0 & 
                                                gse18092.r5.damid.tss$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in\nconservative non-domains", outline = F)
  text(1.5, 2000, paste(length(gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 0 &
                                                                     gse18092.r5.damid.tss$LAM == 0 &
                                                                     gse18092.r5.damid.tss$PC == 0 & 
                                                                     gse18092.r5.damid.tss$truediff == 1]),
                       "genes\np =",
                       format(wilcox.test(gse18092.r5.damid.tss$WT.exp.avg[gse18092.r5.damid.tss$HP1 == 0 &
                                                                             gse18092.r5.damid.tss$LAM == 0 &
                                                                             gse18092.r5.damid.tss$PC == 0 & 
                                                                             gse18092.r5.damid.tss$truediff == 1],
                                          gse18092.r5.damid.tss$HP1KD.exp.avg[gse18092.r5.damid.tss$HP1 == 0 &
                                                                                gse18092.r5.damid.tss$LAM == 0 &
                                                                                gse18092.r5.damid.tss$PC == 0 & 
                                                                                gse18092.r5.damid.tss$truediff == 1],
                                          alt = "l")$p.value, digits = 2)))
  par(op)
dev.off()

boxplot(gse18092.r5.damid$WT.exp.avg[gse18092.r5.damid$truediff == 1],
        gse18092.r5.damid$HP1KD.exp.avg[gse18092.r5.damid$truediff == 1],
        names = c("WT", "HP1 KD"))

wilcox.test(gse18092.r5.damid$WT.exp.avg[gse18092.r5.damid$truediff == 1],
             gse18092.r5.damid$HP1KD.exp.avg[gse18092.r5.damid$truediff == 1], alt = "l")

# 2. Whole Genes

# 2.1. All genes

pdf("all.whole.genes.in.domains.wt.hp1kd.pdf", width = 12)
  op <- par(mfrow = c(1, 4))
  # HP1
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in HP1 domains", outline = F)
  text(1.5, 12, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1]),
                      "genes\np =",
                      format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1],
                                         gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1],
                                         alt = "l")$p.value, digits = 2)))
  # LAM
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in LAM domains", outline = F)
  text(1.5, 11.7, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1]),
                     "genes\np =",
                     format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1],
                                        gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1],
                                        alt = "l")$p.value, digits = 2)))
  # PC
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$PC == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in PC domains", outline = F)
  text(1.5, 12, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1]),
                       "genes\np =",
                       format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$PC == 1],
                                          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1],
                                          alt = "l")$p.value, digits = 2)))
  # Conservative non-domains (WRONG, see diff genes)
  # boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$nondom == 1],
  #         gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$nondom == 1],
  #         names = c("WT", "HP1 KD"),
  #         main = "Gene expression in\nconservative non-domains", outline = F)
  # text(1.5, 9, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$nondom == 1]),
  #                       "genes\np =",
  #                       format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$nondom == 1],
  #                                          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$nondom == 1],
  #                                          alt = "g")$p.value, digits = 2)))
  par(op)
dev.off()

# 2.2. Diff expressed genes

# unlog expr

gse18092.r5.damid.wh.genes <- gse18092.r5.damid.wh.genes %>% mutate(WT.exp.avg = 2^WT.exp.avg,
                                                                    HP1KD.exp.avg = 2^HP1KD.exp.avg) %>%
  filter(!(ID %in% trash_id))

pdf("difex.whole.genes.in.domains.wt.hp1kd.unlog.pdf", width = 12)
  op <- par(mfrow = c(1, 4))
  # HP1
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in HP1 domains", outline = F)
  text(1.5, 2000, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1]),
                      "genes\np =",
                      format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                         gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                         alt = "l")$p.value, digits = 2)))
  # LAM
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in LAM domains", outline = F)
  text(1.5, 1600, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1 & gse18092.r5.damid.wh.genes$truediff == 1]),
                       "genes\np =",
                       format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$LAM == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                          alt = "l")$p.value, digits = 2)))
  # PC
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$PC == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in PC domains", outline = F)
  text(1.5, 1500, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1 & gse18092.r5.damid.wh.genes$truediff == 1]),
                       "genes\np =",
                       format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$PC == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$PC == 1 & gse18092.r5.damid.wh.genes$truediff == 1],
                                          alt = "l")$p.value, digits = 2)))
  # Conservative non-domains
  boxplot(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 0 &
                                                  gse18092.r5.damid.wh.genes$LAM == 0 &
                                                  gse18092.r5.damid.wh.genes$PC == 0 &
                                                  gse18092.r5.damid.wh.genes$truediff == 1],
          gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 0 &
                                                     gse18092.r5.damid.wh.genes$LAM == 0 &
                                                     gse18092.r5.damid.wh.genes$PC == 0 &
                                                     gse18092.r5.damid.wh.genes$truediff == 1],
          names = c("WT", "HP1 KD"),
          main = "Gene expression in\nconservative non-domains", outline = F)
  text(1.5, 3000, paste(length(gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 0 &
                                                                        gse18092.r5.damid.wh.genes$LAM == 0 &
                                                                        gse18092.r5.damid.wh.genes$PC == 0 &
                                                                        gse18092.r5.damid.wh.genes$truediff == 1]),
                        "genes\np =",
                        format(wilcox.test(gse18092.r5.damid.wh.genes$WT.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 0 &
                                                                                   gse18092.r5.damid.wh.genes$LAM == 0 &
                                                                                   gse18092.r5.damid.wh.genes$PC == 0 &
                                                                                   gse18092.r5.damid.wh.genes$truediff == 1],
                                           gse18092.r5.damid.wh.genes$HP1KD.exp.avg[gse18092.r5.damid.wh.genes$HP1 == 0 &
                                                                                      gse18092.r5.damid.wh.genes$LAM == 0 &
                                                                                      gse18092.r5.damid.wh.genes$PC == 0 &
                                                                                      gse18092.r5.damid.wh.genes$truediff == 1],
                                           alt = "l")$p.value, digits = 2)))
  par(op)
dev.off()


# Levels of expression of dif expressed (according to expression MA) genes (whole bodies)
# residing in HP1 domains according to RNA-seq in Kc
dif.genes <- gse18092.r5.damid.wh.genes$CLONE_ID_LIST[gse18092.r5.damid.wh.genes$truediff == 1]

dif.genes.in.hp1 <- gse18092.r5.damid.wh.genes$CLONE_ID_LIST[gse18092.r5.damid.wh.genes$HP1 == 1 & gse18092.r5.damid.wh.genes$truediff == 1]

rna.seq.data <- new.env()
load("../RNAseq_vs_DamID/RNAseq_DamID.RData", envir = rna.seq.data)

expr.data.kc <- rna.seq.data$kc.rsem %>% filter(id %in% dif.genes)

expr.data.kc$in.hp1.ma.damid <- ifelse(expr.data.kc$id %in% dif.genes.in.hp1, "bound", "not bound")
pdf("kc.dif.expressed.genes.HP1.bound.RNA-seq.pdf", width = 4)
  boxplot(TPM ~ in.hp1.ma.damid, data = expr.data.kc, outline = F,
        ylab = "TPM",
        main = "Expression of HP1 targets\nfrom differentially expressed\ngenes in KC")
dev.off()

# boxplot(expr.data.kc$TPM, outline = F)

wilcox.test(expr.data.kc$TPM[expr.data.kc$in.hp1.ma.damid == "bound"],
            expr.data.kc$TPM[expr.data.kc$in.hp1.ma.damid != "bound"],
            alt = "l")

expr.data.kc %>% group_by(in.hp1.ma.damid) %>% summarise(n = n())

# The same but for TSS

dif.tss <- gse18092.r5.damid.tss$CLONE_ID_LIST[gse18092.r5.damid.tss$truediff == 1]
dif.tss.in.hp1 <- gse18092.r5.damid.tss$CLONE_ID_LIST[gse18092.r5.damid.tss$HP1 == 1 & gse18092.r5.damid.tss$truediff == 1]

expr.data.kc <- rna.seq.data$kc.rsem %>% filter(id %in% dif.tss)

expr.data.kc$in.hp1.ma.damid <- ifelse(expr.data.kc$id %in% dif.tss.in.hp1, "bound", "not bound")
pdf("kc.dif.expressed.tss.HP1.bound.RNA-seq.pdf", width = 4)
boxplot(TPM ~ in.hp1.ma.damid, data = expr.data.kc, outline = F,
        ylab = "TPM",
        main = "Expression of HP1 targets\nfrom differentially expressed\ngenes in KC")
dev.off()

wilcox.test(expr.data.kc$TPM[expr.data.kc$in.hp1.ma.damid == "bound"],
            expr.data.kc$TPM[expr.data.kc$in.hp1.ma.damid != "bound"],
            alt = "l")


expr.data.kc %>% group_by(in.hp1.ma.damid) %>% summarise(n = n())

gse18092.damid.diff.ex.wh.genes <- gse18092.r5.damid.wh.genes[gse18092.r5.damid.wh.genes$truediff == 1,] %>% 
  dplyr::select(-truediff, -TSS, -Sequence.Source, -Target.Description)

# remove fbgn dups
dups <- unique(gse18092.damid.diff.ex.wh.genes$CLONE_ID_LIST[duplicated(gse18092.damid.diff.ex.wh.genes$CLONE_ID_LIST)])
dups.df <- gse18092.damid.diff.ex.wh.genes[gse18092.damid.diff.ex.wh.genes$CLONE_ID_LIST %in% dups,] %>% arrange(CLONE_ID_LIST)
trash_id <- dups.df$ID[-c(1,5,6,9,10,13,14,17,18,20,22,24,26)]

gse18092.damid.diff.ex.wh.genes <-  gse18092.damid.diff.ex.wh.genes %>% filter(!(ID %in% trash_id) )

write.table(gse18092.damid.diff.ex.wh.genes, "diff.genes.whole.domains.euc.csv", quote = F, dec = ",",
            sep = ";", row.names = F)

# save.image("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/HP1.kd.MA.Kc/HP1.kd.ma.RData")

# Compare rna-seq and rna ma

rsem.ma.wt <- merge(expr.data.kc, gse18092.r5.damid.wh.genes[,c(6,12,16)], by.x = "id", by.y = "CLONE_ID_LIST")


shev.data <- fread("final microarray data joint with RSEM_in domains.csv")
gse18092.r5.twof <- gse18092.r5[gse18092.r5$twofold == 1,]
gse18092.r5.twof.dif <- gse18092.r5[gse18092.r5$twofold == 1 & gse18092.r5$truediff == 1,]
gse18092.r5.twof.ndif <- gse18092.r5[gse18092.r5$twofold == 1 & gse18092.r5$truediff == 0,]

