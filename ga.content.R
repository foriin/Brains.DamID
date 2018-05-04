library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

rm(list=ls())

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/ga.cont/")

euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

###################################
# DamID data
damid.fol <- "~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/BioHMM2.qn.full.PC.HMM3/"
domains <- new.env()
for (i in dir(damid.fol)[grepl("LAM", dir(damid.fol))]){
  df <- fread(file.path(damid.fol, i), skip = 1, col.names = c("chr", "start", "end")) %>% filter(chr %in% euc)
  assign(sub("\\.domains\\.bed", "", i),
         GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
         envir = domains)
}
rm(i)

lams <- list()

for (i in 1:5){
  lams[i] <- get(ls(envir = domains)[i], envir = domains)
}

cLADs <- Reduce(GenomicRanges::intersect, lams)

cLADs.df <- data.frame(
  chr = seqnames(cLADs),
  start = start(cLADs),
  end = end(cLADs)
) %>% mutate(chr = sub("chr", "", chr))

cat(paste(paste(cLADs.df$chr, cLADs.df$start, sep = ":"), cLADs.df$end, sep = ".."), file = "cLADs.coords.txt", sep = "\n")

# Now go at http://flybase.org/static_pages/downloads/COORD.html and convert this coordinates to r6 release
# because batch download (http://flybase.org/static_pages/downloads/ID.html) works only with r6 release
# then retrieve fasta for those coordinates

ciLADs <- Reduce(GenomicRanges::intersect, lapply(lams, gaps))

ciLADs.df <- data.frame(
  chr = seqnames(ciLADs),
  start = start(ciLADs),
  end = end(ciLADs)
) %>% mutate(chr = sub("chr", "", chr))

cat(paste(paste(ciLADs.df$chr, ciLADs.df$start, sep = ":"), ciLADs.df$end, sep = ".."), file = "ciLADs.coords.txt", sep = "\n")

detach("package:dplyr", unload=TRUE)
fLADs <- setdiff(Reduce(union, lams), cLADs)
library(dplyr)

fLADs.df <- data.frame(
  chr = seqnames(fLADs),
  start = start(fLADs),
  end = end(fLADs)
) %>% mutate(chr = sub("chr", "", chr))

cat(paste(paste(fLADs.df$chr, fLADs.df$start, sep = ":"), fLADs.df$end, sep = ".."), file = "fLADs.coords.txt", sep = "\n")


Kc.df <- data.frame(
  chr = seqnames(domains$Kc167.LAM),
  start = start(domains$Kc167.LAM),
  end = end(domains$Kc167.LAM)
) %>% mutate(chr = sub("chr", "", chr))

cat(paste(paste(Kc.df$chr, Kc.df$start, sep = ":"), Kc.df$end, sep = ".."), file = "Kc.coords.txt", sep = "\n")

Kc.gaps <- gaps(domains$Kc167.LAM)

Kc.gaps.df <- data.frame(
  chr = seqnames(Kc.gaps),
  start = start(Kc.gaps),
  end = end(Kc.gaps)
) %>% mutate(chr = sub("chr", "", chr))

cat(paste(paste(Kc.gaps.df$chr, Kc.gaps.df$start, sep = ":"), Kc.gaps.df$end, sep = ".."), file = "Kc.gaps.coords.txt", sep = "\n")
