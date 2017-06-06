library("dplyr")
library("data.table")
library(GenomicRanges)
rm(list = ls())
setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/02.04.17_all_samples_300/CSV/")
dam.norm <- fread("07.Normalized.MeanDATA.csv")

hmms <- list.files("../BioHMM/")

e1 <- new.env(parent = baseenv())

for (i in hmms){
  y <- sub("(^[^.]+\\.[^.]+).*", "\\1", i)
  ynas <- cbind(dam.norm[, 2:4], subset(dam.norm, select = grep(y, names(dam.norm))))
  ynas <- ynas[!complete.cases(ynas),]
  names(ynas)[ncol(ynas)] <- "score"
  # ynas$chr <- paste0("chr", ynas$chr)
  ycomp <- fread(paste0("../BioHMM/", i), skip = 1)
  names(ycomp) <- c("chr", "start", "end")
  ycomp$score <- 1
  yall <- rbind(ycomp, ynas) %>% arrange(chr, start) %>% mutate(add = 0)
  for (i in 2:(nrow(yall) - 1)){
    yall$add[i] <- yall$score[i - 1] + yall$score[i + 1] + yall$start[i] - yall$end[i-1] + 
      yall$start[i + 1] - yall$end[i]
  }
  yall <- yall %>% filter(!is.na(score) | add < 4)
  assign(y, yall, envir = e1)
}

dir.create("../hmm3")
for (bedranges in ls(envir = e1)){
  bed <- get(bedranges, envir = e1)
  gr <- GRanges(seqnames = Rle(bed$chr),  ranges = IRanges(start = bed$start, end = bed$end)) %>% 
    reduce() 
  write.table(data.frame(chr = seqnames(gr), start = start(gr), end = end(gr)),
              paste0("../hmm3/", bedranges, ".bed"), row.names = F, quote = F, col.names = F,
              sep = '\t')
  
}

