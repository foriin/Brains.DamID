# Libs
library(dplyr)
library(GenomicRanges)
library(data.table)
library(regioneR)

# Используем системный вызов bedtools-shuffle для рандомной перестановки доменов,
# т.к. все что удалось найти в биокондукторе и т.п. экскрементивного качества

rm(list=ls())

# Functions
bedTools.shuffle <- function(bed.gr, opt.string="-chrom"){
  bed <- data.frame(chr = seqnames(bed.gr),
                    start = start(bed.gr),
                    end = end(bed.gr))
  
  bed.file <- tempfile()
  out <- tempfile()
  options(scipen = 99)
  
  write.table(bed, file = bed.file, quote = F, sep = "\t", col.names = F, row.names = F)
  command = paste("bedtools shuffle -i", bed.file,
                  "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, ">", out)
  # cat(command, "\n")
  try(system(command))
  
  res=read.table(out, header = F)
  colnames(res) <- c("chr", "start", "end")
  unlink(bed.file); unlink(out)
  return(GRanges(seqnames = Rle(res$chr), ranges = IRanges(start = res$start, end = res$end)))
}

################################



setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/Ubiq_genes/")

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

ub.genes <- read.csv("ubiqex.genes.csv", stringsAsFactors = F)
names(ub.genes)[4] <- "chr"
ub.genes$chr <- paste0("chr", ub.genes$chr)
ub.genes <- ub.genes %>% filter(chr %in% euc)

ub.gr <- GRanges(seqnames = Rle(ub.genes$chr),
                 ranges = IRanges(start = ub.genes$promoter.position, width = 1))
seqlengths(ub.gr) <- chr.len[names(seqlengths(ub.gr))]

# see how many genes on each chromosome
# lapply(split(ub.gr, seqnames(ub.gr)), length)

dir.w.beds <- paste0("~/IMG/Projects/",
                     "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
                     "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/BioHMM2.qn.full.PC.HMM3/")

beds <- new.env()
for (i in dir(dir.w.beds)[grepl("bed", dir(dir.w.beds))]){
  assign(sub("\\.domains.bed", "", i), fread(file.path(dir.w.beds, i), skip = 1), envir = beds)
}
rm(i)

prots <- unique(sub("[^.]+\\.([^.]+)", "\\1", names(beds)))

doms.adoms <- lapply(prots, function(prot){
  retl <- list()
  beds.p.i <- which(grepl(prot,  names(beds)))
  prots.grs <- lapply(beds.p.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr %in% euc)
    GRanges(seqnames = Rle(df$chr),
            ranges = IRanges(start = df$start, end = df$end))
  })
  names(prots.grs) <- sub("\\..+$", "", names(beds)[beds.p.i])
  prots.grs <- prots.grs[order(names(prots.grs))]
  
  # 5x intersect
  tot.intrsct <- Reduce(GenomicRanges::intersect, prots.grs)
  seqlengths(tot.intrsct) <- chr.len[names(seqlengths(tot.intrsct))]
  retl[["cons.doms"]] <- tot.intrsct
  
  # anti-domains
  prots.gaps <- lapply(prots.grs, gaps)
  gaps.intrsct <- Reduce(GenomicRanges::intersect, prots.gaps)
  seqlengths(gaps.intrsct) <- chr.len[names(seqlengths(gaps.intrsct))]
  retl[["cons.antidoms"]] <- gaps.intrsct
  return(retl)
})

names(doms.adoms) <- prots



numOverlaps(doms.adoms[[1]]$cons.antidoms, ub.gr)
x <- randomizeRegions(doms.adoms[[1]]$cons.doms, genome = "dm3", allow.overlaps = F)
sum(width(reduce(doms.adoms[[1]]$cons.doms)))
sum(width(reduce(x)))
seqlengths(x)


domains.p.t.all.ran <- sapply(prots, function(pr) {
  print(pr)
  # Genes in conservative domains
  ub.genes.in.cons.dom <- GenomicRanges::intersect(doms.adoms[[pr]]$cons.doms, ub.gr)
  # dom.genes <- merge(ub.genes, data.frame(chr = seqnames(ub.genes.in.cons.dom),
  #                                         promoter.position = start(ub.genes.in.cons.dom)),
  #                    by = c("chr", "promoter.position"))$CG
  in.dom.ratio <- length(ub.genes.in.cons.dom)/length(ub.gr)
  dom.perm.stat <- sapply(1:10000, function(i){
    length(GenomicRanges::intersect(bedTools.shuffle(doms.adoms[[pr]]$cons.doms),
                                    bedTools.shuffle(ub.gr))) / length(ub.gr)
    
  })
  pval.d <- sum(dom.perm.stat <= in.dom.ratio)/10000
  print(pval.d)
  
  ub.genes.in.cons.antidom <- GenomicRanges::intersect(doms.adoms[[pr]]$cons.antidoms, ub.gr)
  # adom.genes <- merge(ub.genes, data.frame(chr = seqnames(ub.genes.in.cons.antidom),
  #                                         promoter.position = start(ub.genes.in.cons.antidom)),
  #                    by = c("chr", "promoter.position"))$CG
  in.adom.ratio <- length(ub.genes.in.cons.antidom)/length(ub.gr)
  adom.perm.stat <- sapply(1:10000, function(i){
    length(GenomicRanges::intersect(bedTools.shuffle(doms.adoms[[pr]]$cons.antidoms),
                                    bedTools.shuffle(ub.gr))) / length(ub.gr)
    
  })
  pval.a <- sum(adom.perm.stat >= in.adom.ratio)/10000
  print(pval.a)
  c("In cons domains" = in.dom.ratio, "dom p-value" = pval.d,
    "In cons adomains" = in.adom.ratio, "adom p-vlaue" = pval.a)
  # ub.genes[[pr]] <- NA
  # ub.genes[ub.genes$CG %in% dom.genes,][[pr]] <- 1
  # ub.genes[ub.genes$CG %in% adom.genes,][[pr]] <- 0
  # 
  # 
 
  # 
  # 
  # cat(paste0(pr,
  #            "\nin domains: ", length(dom.genes)/nrow(ub.genes)*100, "%\n",
  #            "in anti-domains: ", round(length(adom.genes)/nrow(ub.genes)*100, digits = 2), "%\n"),
  #            file = "ubiquitous.expr.genes.in.domains.txt",
  #            append = T)
})


x <- bedTools.shuffle(doms.adoms[[1]]$cons.doms)

sum(width(reduce(x)))

sum(width(reduce(doms.adoms[[1]]$cons.doms)))
sum(width(reduce(shuffle(doms.adoms[[1]]$cons.doms))))


write.table(ub.genes, "ub.genes.in.cons.domains.csv", quote = F, sep = ";", row.names = F)

write.table(as.data.frame(t(domains.p.t.all.ran)), "ubiq.genes.in.domains.pvalues.csv", quote = F, sep = ";", row.names = F)


lapply(doms.adoms, function(lst){
  sapply(lst, seqlengths)
})


