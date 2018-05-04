library(xlsx)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(reshape2)

setwd("/home/artem/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/neuro_and_glia_RNAseq/")


bg.dir <- paste0("~/IMG/Projects/",
                      "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
                      "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/Bedgraph/")
csv.dir <- paste0("~/IMG/Projects/",
                 "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
                 "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/CSV/")

damid.df <- fread(file.path(csv.dir, "07.Normalized.MeanDATA.csv"))
damid.df[, index := seq_len(.N), by = chr]
names(damid.df) <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", names(damid.df))
damid.df <- damid.df %>% select(-grep("FB|BR", names(damid.df))) %>% mutate(chr = paste0("chr", chr))

bg <- new.env()

for (i in names(damid.df)[5:(ncol(damid.df) - 1)]){
  assign(i,
         makeGRangesFromDataFrame(
           data.frame(damid.df %>% select(c(2:4, ncol(damid.df))), score = damid.df[[i]]),
           keep.extra.columns = T),
         envir = bg)
}

x <- findOverlaps(glia.tss.gr, bg$Glia.LAM)

lol <- sapply(1:length(x), function(ind){
  tss.n <- x@from[ind]
  tss.c <- x@to[ind]
  if (glirnaseq[tss.n, 7] == 1){
    bg$LAM.Glia$score[c(tss.c - 1, tss.c, tss.c + 1)]
  } else {
    bg$LAM.Glia$score[c(tss.c + 1, tss.c, tss.c - 1)]
  }
})

apply(lol, 1, median)
median(lol[1,])

damid.tss <- function(rnaseq, damid){
  # rnaseq <- data.frame with tissue-specific
  # genes, with TPMs cut in three categories
  # via `cut(rnaseq$TPM, quantile(rnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
  # labels = c("weak", "mid", "high"), include.lowest = T)` command
  # strand information also necessary

  sapply(levels(rnaseq$TPM.cut), function(ct){
    dt.ct <- rnaseq %>% filter(TPM.cut == ct)
    ct.tss.gr <- GRanges(
      seqnames = Rle(dt.ct$chr),
      ranges = IRanges(
        start = as.integer(dt.ct$start),
        width = 1
      )
    )
    overlaps <- findOverlaps(ct.tss.gr, damid)
    vic <- sapply(1:length(overlaps), function(ind){
      tss.n <- overlaps@from[ind]
      tss.c <- overlaps@to[ind]
      if (dt.ct[tss.n, 7] == 1){
        damid$score[-2:2 + tss.c]
      } else {
        damid$score[2:-2 + tss.c]
      }
    })
    # print(head(vic))
    apply(vic, 1, median)
  })
  # print(head(as_data_frame(vic)))
  
  # apply(vic, 1, median)
  
}


# sapply(a, function(x) apply(x, 1,  median))
a <- damid.tss(glirnaseq, bg$Glia.HP1) %>% as.data.frame()
a$tss <- -2:2
b <- melt(a,measure.vars = 1:3)

ggplot(b, aes(x = tss, y = value))+
  geom_line(aes(col = variable))+
  ylim(-2, 1)



damid.tss(glirnaseq, bg$HP1.Glia) %>% as.data.frame()

damid.tss(glirnaseq, bg$PC.Glia) %>% as.data.frame() 

damid.tss(neurnaseq, bg$LAM.NRN) %>% as.data.frame() 






damid.tss.2 <- function(rnaseq, damid, range = 2, fun = median, ...){
  # rnaseq <- data.frame with tissue-specific
  # genes, with TPMs cut in three categories
  # via `cut(rnaseq$TPM, quantile(rnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
  # labels = c("weak", "mid", "high"), include.lowest = T)` command
  # strand information also necessary
  
  sapply(levels(rnaseq$TPM.cut), function(ct){
    dt.ct <- rnaseq %>% filter(TPM.cut == ct)
    ct.tss.gr <- GRanges(
      seqnames = Rle(dt.ct$chr),
      ranges = IRanges(
        start = as.integer(dt.ct$TSS),
        width = 1
      )
    )
  
    overlaps <- findOverlaps(ct.tss.gr, damid)
    print(length(overlaps))
    # print(length(which(damid[overlaps@to]$index < 10)))
    vic <- sapply(1:length(overlaps@to), function(ind){
      tss.n <- overlaps@from[ind]
      tss.c <- overlaps@to[ind]
      if (dt.ct[tss.n, 7] == 1){
        damid$score[-range:range + tss.c]
      } else {
        damid$score[range:-range + tss.c]
      }
    })

    # print(head(vic))
    if (range > 0){
      return(apply(vic, 1, fun, ...))
    }else{
      return(fun(vic, ...))
    }
  })
   # print(head(as_data_frame(vic)))
  
}

xx <- damid.tss.2(glirnaseq, bg$Glia.HP1, range = 5,na.rm = T)

# Make plots!

prots <- c("LAM", "HP1", "PC")

for (pr in prots){
  rng <- 5
  a <- damid.tss.2(glirnaseq %>% filter(HP1 == 1), bg[[paste0("Glia.", pr)]], range = rng, na.rm =T) %>% as.data.frame()
  a$tss <- factor(-rng:rng)
  b <- melt(a,measure.vars = 1:3)
  p <- ggplot(b, aes(x = tss, y = value))+
    geom_line(aes(col = variable, group = variable))+
    ylim(-2, 1)+
    theme_bw()
  pdf(paste0(pr, ".Glia.tss.dist.hp1.pdf"))
    print(p)
  dev.off()
}

for (pr in prots){
  rng <- 5
  a <- damid.tss.2(neurnaseq %>% filter(HP1 == 1), bg[[paste0("NRN.", pr)]], range = rng, na.rm = T) %>% as.data.frame()
  a$tss <- factor(-rng:rng)
  b <- melt(a,measure.vars = 1:3)
  p <- ggplot(b, aes(x = tss, y = value))+
    geom_line(aes(col = variable, group = variable))+
    ylim(-2, 1)+
    theme_bw()
  pdf(paste0(pr, ".NRN.tss.dist.hp1.pdf"))
  print(p)
  dev.off()
}
 
# Save to excel!

tss.gli <- createWorkbook()

sheet1 <- createSheet(tss.gli, "LAM")
sheet2 <- createSheet(tss.gli, "HP1")
sheet3 <- createSheet(tss.gli, "PC")

addDataFrame(damid.tss.2(glirnaseq, bg$Glia.LAM) %>% as.data.frame() %>%
               mutate(tss = -5:5), row.names = F,
             sheet1)
addDataFrame(damid.tss.2(glirnaseq, bg$Glia.HP1) %>% as.data.frame() %>%
               mutate(tss = -5:5), row.names = F,
             sheet2)
addDataFrame(damid.tss.2(glirnaseq, bg$Glia.PC) %>% as.data.frame() %>%
               mutate(tss = -5:5), row.names = F,
             sheet3)
saveWorkbook(tss.gli, "Glia_damid_enrichment_in_tss_vicinities.xlsx")


tss.neu <- createWorkbook()

sheet1 <- createSheet(tss.neu, "LAM")
sheet2 <- createSheet(tss.neu, "HP1")
sheet3 <- createSheet(tss.neu, "PC")

addDataFrame(damid.tss.2(neurnaseq, bg$NRN.LAM) %>% as.data.frame() %>% 
               mutate(tss = -5:5), row.names = F,
             sheet1)
addDataFrame(damid.tss.2(neurnaseq, bg$NRN.HP1) %>% as.data.frame() %>%
               mutate(tss = -5:5), row.names = F,
             sheet2)
addDataFrame(damid.tss.2(neurnaseq, bg$NRN.PC) %>% as.data.frame() %>%
               mutate(tss = -5:5), row.names = F,
             sheet3)

saveWorkbook(tss.neu, "Neurons_damid_enrichment_in_tss_vicinities.xlsx")
