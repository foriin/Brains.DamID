library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(GenomicRanges)
library(tximport)

rm(list=ls())

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/")
genes <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf") %>% 
  select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('.*gene_id "(FBgn[0-9]+)";.*', '\\1', attr), gene_name = sub('.*gene_name "([^;]+)";.*', '\\1', attr),
         tss = ifelse(strand == "+", start, end)) %>% 
  select(-attr)
genes$chr <- paste0("chr", genes$chr)

rna.dat <- merge(read.table("salmon.by.gene.csv", header = T), genes, by = "id")
names(rna.dat)[2:3] <- c("brains.exp", "fatbody.exp")
###################################
# DamID data
setwd("../DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/02.04.17_all_samples_300/BioHMM2/intersected/")
domains <- new.env()
for (i in dir()[grepl("BR|FB", dir())]){
   df <- fread(i, skip = 1, col.names = c("chr", "start", "end")) 
   assign(sub("\\.domains\\.intersected\\.bed", "", i),
          GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
          envir = domains)
}
rm(i)

tss.gr <- GRanges(seqnames = Rle(rna.dat$chr),  ranges = IRanges(start = rna.dat$tss, width = 1, names = rna.dat$id))

rna.dat<- cbind(rna.dat, sapply(ls(envir = domains), function(grname){
  gr <- get(grname, envir = domains)
  in_domain <- subsetByOverlaps(tss.gr, gr)@ranges@NAMES
  ifelse(rna.dat$id %in% in_domain, 1, 0)
}))



boxplot(rna.dat$brains.exp[rna.dat$BR.PC == 1 & rna.dat$brains > 0], rna.dat$brains.exp[rna.dat$BR.PC == 0 & rna.dat$brains.exp > 0], outline = F)

boxplot(rna.dat$fatbody.exp[rna.dat$FB.HP1 == 1 & rna.dat$fatbody.exp > 0], rna.dat$fatbody.exp[rna.dat$FB.HP1 == 0 & rna.dat$fatbody.exp > 0], outline = F)

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/")
write.table(rna.dat[,c(1,8,4,5,6,9,7,2,3,10:15)], "BR_FB_rnaseq_vs_damid.csv",
            row.names = F, quote = F, col.names = T, dec = ",",
            sep = '\t')

