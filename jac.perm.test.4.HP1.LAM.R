# Libs
library(dplyr)
library(GenomicRanges)
library(data.table)

rm(list = ls())

bedTools.shuffle.jac <- function(bed.1, bed.2, shuf = F, opt.string="-chrom"){
  
  bed.file.1 <- tempfile()
  bed.file.2 <- tempfile()
  
  shuf.1 <- tempfile()
  shuf.2 <- tempfile()
  
  jac <- tempfile()
  
  options(scipen = 99)
  
  write.table(bed.1, file = bed.file.1, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(bed.2, file = bed.file.2, quote = F, sep = "\t", col.names = F, row.names = F)
  if (shuf){
    command = paste("bedtools shuffle -i", bed.file.1,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.1, ";",
                    "bedtools shuffle -i", bed.file.2,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.2, ";",
                    "bedtools jaccard -a", shuf.1, "-b", shuf.2, ">", jac)
    # cat(command, "\n")
  } else {
    command = paste("bedtools jaccard -a", bed.file.1, "-b", bed.file.2, ">", jac)
  }
  try(system(command))
  
  res=read.table(jac, header = T)
  unlink(bed.file.1); unlink(bed.file.2); unlink(shuf.2); unlink(shuf.1); unlink(jac)
  return(res$jaccard)
}

setwd(paste0("~/IMG/Projects/",
       "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
       "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/BioHMM2.qn.full.PC.HMM3/"))

# load("jac.perm.test.hp1.lam.RData")

chr.len <- c("chr2L" = 23011544,
             "chr2R" = 21146708,
             "chr3L" = 24543557,
             "chr3R" = 27905053,
             "chr4" = 1351857,
             "chrX" = 22422827,
             "chr2LHet" = 368872,
             "chr2RHet" = 3288761,
             "chr3LHet" = 2555491,
             "chr3RHet" = 2517507,
             "chrXHet" = 204112,
             "chrYHet" = 347038)

euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)

prots <- unique(sub("[^.]+\\.([^.]+)", "\\1", names(beds)))
prots.n <- prots[-2]

tiss <- unique(sub("^([^.]+)\\.([^.]+)", "\\1", names(beds)))

jac.perm.pv.less <- sapply(tiss[-c(1,2,3)], function(tis){
  tiss.hp1.lam.doms <- names(beds)[grepl(tis, names(beds)) & grepl(paste(prots.n, collapse = "|"), names(beds))]
  jac <- bedTools.shuffle.jac(get(tiss.hp1.lam.doms[1], envir = beds),
                              get(tiss.hp1.lam.doms[2], envir = beds))
  print(jac)
  jac.shuf <- sapply(1:10000, function(i){
    bedTools.shuffle.jac(get(tiss.hp1.lam.doms[1], envir = beds),
                         get(tiss.hp1.lam.doms[2], envir = beds), shuf = T)
  })
  
  pval <- sum(jac.shuf <= jac)/10000
  print(pval)
  return(c(tis, jac, pval))
  
  })

jac.out <- rbind(as.data.frame(t(jac.perm.pv), stringsAsFactors = F),
                 as.data.frame(t(jac.perm.pv.less), stringsAsFactors = F)[1,]) 
names(jac.out) <- c("Tissue", "Jaccard", "pval")
jac.out$Jaccard <- as.numeric(jac.out$Jaccard)
jac.out$pval <- as.numeric((jac.out$pval))

write.table(jac.out, "perm.test.results.for.hp1.lam.n.10000.csv", quote = F, dec = ",", sep = "\t", row.names = F)
