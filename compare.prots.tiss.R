library("dplyr")
library("data.table")
library(GenomicRanges)
rm(list = ls())

setwd(paste0("~/IMG/Projects/",
             "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/21.07.17_1kbins_NRN/",
             "sum_rep_qn/BioHMM2.qn.full.PC.HMM3/"))

frac.domains.euc <- function(df1, df2){
  df1 <- df1 %>% filter(chr %in% c("2L", "2R", "3L", "3R", "X"))
  df2 <- df2 %>%  filter(chr %in% c("2L", "2R", "3L", "3R", "X"))
  round(sum(df2$end - df2$start) / sum(df1$end - df1$start), digits = 3)
}

beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.intersected.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)

prots <- unique(sub("[^.]+\\.([^.]+)\\..*", "\\1", names(beds)))

protcom <- lapply(prots, function(pr){
  prot.gr <- lapply(names(beds)[grepl(pr, names(beds))], function(nom){
    bed <- get(nom, envir = beds)
    names(bed) <- c("chr", "start", "end")
    gr <- GRanges(seqnames = Rle(bed$chr),
                  ranges = IRanges(start = bed$start, end = bed$end))
    })
  intersect.gr <- Reduce(intersect, prot.gr)
  print(names(prot.gr))
  data.frame(chr = seqnames(intersect.gr), start = start(intersect.gr), end = end(intersect.gr))

})

names(protcom) <- prots

lapply(names(protcom), function(nim){
  txt.name <- paste0(nim, ".domains.proportion.txt")
  lapply(names(beds)[grepl(nim, names(beds))], function(nom){
    bed <- get(nom, envir = beds)
    names(bed) <- c("chr", "start", "end")
    print(frac.domains(bed, protcom[[nim]]))
    cat(sub("bed", "", gsub("\\.", " ", nom)), "\t", frac.domains(bed, protcom[[nim]]), "\n",
        file = file.path("./ProtCom", txt.name), append = T)
  })
})

dir.create("ProtCom", showWarnings = F)
lapply(prots, function(prot){
  write.table(protcom[[prot]], file=file.path("ProtCom", paste0(prot, ".total.intersect.bed")),
                                              sep="\t", row.names=F, col.names=F, quote=F,
                                              dec=".", append=T)
})

tiss <- unique(sub("(^[^.]+)\\..*", "\\1", names(beds)))

tisscom <- lapply(tiss, function(tis){
  prot.gr <- lapply(names(beds)[grepl(tis, names(beds))], function(nom){
    bed <- get(nom, envir = beds)
    names(bed) <- c("chr", "start", "end")
    gr <- GRanges(seqnames = Rle(bed$chr),
                  ranges = IRanges(start = bed$start, end = bed$end))
  })
  intersect.gr <- Reduce(intersect, prot.gr)
  tot.tiss.in <- data.frame(chr = seqnames(intersect.gr), start = start(intersect.gr), end = end(intersect.gr))
  
})

names(tisscom) <- tiss

dir.create("TissCom", showWarnings = F)
lapply(tiss, function(tis){
  write.table(tisscom[[tis]], file=file.path("TissCom", paste0(tis, ".total.intersect.bed")),
              sep="\t", row.names=F, col.names=F, quote=F,
              dec=".", append=T)
})

lapply(names(tisscom), function(nim){
  txt.name <- paste0(nim, ".domains.proportion.txt")
  lapply(names(beds)[grepl(nim, names(beds))], function(nom){
    bed <- get(nom, envir = beds)
    names(bed) <- c("chr", "start", "end")
    cat(sub("bed", "", gsub("\\.", " ", nom)), "\t", frac.domains(bed, tisscom[[nim]]), "\n",
        file = file.path("./TissCom", txt.name), append = T)
  })
})