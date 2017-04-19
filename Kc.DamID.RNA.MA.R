library(data.table)
library(dplyr)
library(tidyr)

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/Kc167_multi_DamID_Filion_2010/")
kc.damid.hmm <- fread("Kc167_HP1.LAM.PC.DamID.HMM.csv")
chrs <- unique(kc.damid.hmm$chromosome)

tss <- fread("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/HP1.kd.MA.Kc/GSE18092.r5.gene.coord.tss.csv")
tss$chr <- paste0("chr", tss$chr)

tss.damid <- lapply(chrs, function(chrom){
  damid.chr <- kc.damid.hmm %>% filter(chromosome == chrom)
  tss.chr <- tss %>% filter(chr == chrom)
  tss.chr$int <- sapply(tss.chr$TSS, function(TSS) findInterval(TSS, damid.chr$start))
  tss.chr$HP1.d <- NA
  tss.chr$Lam.d <- NA
  tss.chr$Pc.d <- NA
  
  for (i in 1:nrow(tss.chr)){
    print(tss.chr$int[i])
    if (tss.chr$int[i]  != 0 && tss.chr[i, 12] < damid.chr[tss.chr[i,13], 4]) {
      tss.chr$HP1.d[i] <- damid.chr[tss.chr[i,]$int, 5]
      tss.chr$Lam.d[i] <- damid.chr[tss.chr[i,]$int, 6]
      tss.chr$Pc.d[i] <- damid.chr[tss.chr[i,]$int, 7]
    }
  }
  return(tss.chr)
})

tss.damid <- do.call("rbind", tss.damid)

tss.damid2 <- tss.damid %>% 
  filter(!is.na(HP1.d)) %>% 
  mutate(HP1.d = factor(HP1.d),
         Lam.d = factor(Lam.d),
         Pc.d = factor(Pc.d))

fbgn.dup <- unique(tss.damid2$FBgn[duplicated(tss.damid2$FBgn)])

thr.out <- sapply(fbgn.dup, function(dup){
  rowin <- which(tss.damid2$FBgn == dup)
  rowin[rowin != which.max(tss.damid2$WT.exp.avg[rowin])]
  })

thr.out <- do.call("c", thr.out)

tss.damid2 <- tss.damid2[-thr.out, ]

write.table(tss.damid2, "tss.in.damid.domains.csv", row.names = F, sep = '\t', dec = ',')

pdf("Kc167.HP1.DamID.RNA.MA.comparison.pdf", width = 8, height = 5)
  boxplot(tss.damid2$WT.exp.avg[tss.damid2$HP1.d == 0],
          tss.damid2$WT.exp.avg[tss.damid2$HP1.d == 1],
          tss.damid2$HP1KD.exp.avg[tss.damid2$HP1.d == 0],
          tss.damid2$HP1KD.exp.avg[tss.damid2$HP1.d == 1],
          names = c("WT expr no HP1", "WT expr HP1", "dsHP1 expr no HP1", "dsHP1 expr HP1"),
          ylab = "average expression", boxwex = 0.5)
dev.off()

pdf("Kc167.Lamin.DamID.RNA.MA.comparison.pdf", height = 5)
  boxplot(tss.damid2$WT.exp.avg[tss.damid2$Lam.d == 0],
          tss.damid2$WT.exp.avg[tss.damid2$Lam.d == 1], outline = F,
          names = c("WT expr no Lam", "WT expr Lam"),
          ylab = "average expression", boxwex = 0.5)
dev.off()


pdf("Kc167.Polycomb.DamID.RNA.MA.comparison.pdf", height = 5)
  boxplot(tss.damid2$WT.exp.avg[tss.damid2$Pc.d == 0],
          tss.damid2$WT.exp.avg[tss.damid2$Pc.d == 1], outline = F,
          names = c("WT expr no Pc", "WT expr Pc"),
          ylab = "average expression", boxwex = 0.5)
dev.off()
