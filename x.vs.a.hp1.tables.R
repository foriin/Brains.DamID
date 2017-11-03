library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(reshape2)
library(BSgenome.Dmelanogaster.UCSC.dm3)

rm(list = ls())


setwd(paste0("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/HP1.bedgraph.plus.Kc"))
euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
genome <- BSgenome.Dmelanogaster.UCSC.dm3
euc.len <- seqlengths(genome)[euc]

# "Boundaries" of euchromatin

euc.gr <- GRanges(
  seqnames = Rle(euc),
  ranges = IRanges(c(1, 1600000, 1, 1, 1), end = c(22000000, 21146708, 22900000, 27905053, 22300000))
)

euc.euc.len <- width(euc.gr)
names(euc.euc.len) <- names(euc.len)

hp1.dom.fol <- paste0("~/IMG/Projects/",
"HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
"DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/",
"final_variant/HP1.bedgraph.plus.Kc/BioHMM.qn.euc/")

domains <- new.env()
# hp1 domains
# for (i in dir(hp1.dom.fol, pattern = "(HP1).*bed")){
#   df <- fread(file.path(hp1.dom.fol, i), skip = 1, col.names = c("chr", "start", "end")) %>% filter(chr %in%euc)
#   assign(sub("\\.domains\\.bed", "", i),
#          intersect(GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
#                    euc.gr),
#          envir = domains)
# }
# rm(i)

# lam domains

dom.dir <- paste0("~/IMG/Projects/",
                  "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
                  "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
                  "BioHMM2.qn.full.PC.HMM3/")

for (i in dir(dom.dir, pattern = "(HP1|LAM).*bed")){
  df <- fread(file.path(dom.dir, i), skip = 1, col.names = c("chr", "start", "end")) %>% filter(chr %in%euc)
  assign(sub("\\.domains\\.bed", "", i),
         intersect(
           GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
           euc.gr),
         envir = domains)
}
rm(i)

# Get domain coverage for X and mean coverage for autosomes
# split separates granges by chromosomes so that you could get total length
# of domains in a given chromosome

cov.hp1 <- sapply(ls(pattern = "HP1", e = domains), function(nym){
  gr <- get(nym, e = domains)
  dom.sp <- split(gr, seqnames(gr))
  cover <- sapply(dom.sp, function(x) sum(width(x))) / euc.euc.len
  round(c("A" = mean(cover[1:4]), cover[5]), digits = 3)
})

in.lad <- sapply(ls(pat = "HP1", e = domains), function(nom){
  lads <- get(ls(envir = domains)[grepl(sub("HP1", "LAM", nom), ls(envir = domains), ignore.case = T)], envir = domains)
  # lads.len <- sapply(split(lads, seqnames(lads)), function(x) sum(width(x)))
  
  gr <- intersect(
    get(nom, envir = domains),
    lads)
  dom.sp <- split(gr, seqnames(gr))
  cover <- sapply(dom.sp, function(x) sum(width(x))) / euc.euc.len
  round(c("A" = mean(cover[1:4]), cover[5]), digits = 3)
  })

out.lad <- sapply(ls(pat = "HP1", e = domains), function(nom){
  out.lads <- gaps(get(ls(envir = domains)[grepl(sub("HP1", "LAM", nom), ls(envir = domains), ignore.case = T)], envir = domains))
  # out.lads.len <- sapply(split(out.lads, seqnames(out.lads)), function(x) sum(width(x)))
  
  gr <- intersect(
    get(nom, envir = domains),
    out.lads)
  dom.sp <- split(gr, seqnames(gr))
  cover <- sapply(dom.sp, function(x) sum(width(x))) / euc.euc.len
  round(c("A" = mean(cover[1:4]), cover[5]), digits = 3)
})

cat("HP1 domains\n", file = "x.vs.a.2.csv")
write.table(cov.hp1, file = "x.vs.a.2.csv", sep = ";", quote = F, dec = ",", row.names = T, col.names = T, append = T)
cat("LADs\n", file = "x.vs.a.2.csv", append = T)
write.table(in.lad, file = "x.vs.a.2.csv", sep = ";", quote = F, dec = ",", row.names = T, col.names = T, append = T)
cat("interLADs\n", file = "x.vs.a.2.csv", append = T)
write.table(out.lad, file = "x.vs.a.2.csv", sep = ";", quote = F, dec = ",", row.names = T, col.names = T, append = T)

# sapply(split(domains$BR.LAM, seqnames(domains$BR.LAM)), function(x) sum(width(x)))



