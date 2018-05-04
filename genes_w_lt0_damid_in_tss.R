

damid.tss.3 <- function(rnaseq, damid, range = 2, fun = median, ...){
  # rnaseq <- data.frame with tissue-specific
  # genes, with TPMs cut in three categories
  # via `cut(rnaseq$TPM, quantile(rnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
  # labels = c("weak", "mid", "high"), include.lowest = T)` command
  # strand information also necessary
  
  tss.gr <- GRanges(
      seqnames = Rle(rnaseq$chr),
      ranges = IRanges(
        start = as.integer(rnaseq$TSS),
        width = 1
      )
    )
    
    overlaps <- findOverlaps(tss.gr, damid)
    print(length(unique(overlaps@from)))
    # print(length(which(damid[overlaps@to]$index < 10)))
    vic <- sapply(1:length(overlaps@to), function(ind){
      tss.n <- overlaps@from[ind]
      tss.c <- overlaps@to[ind]
      if (rnaseq[tss.n, 7] == 1){
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
}

damid.tss.3(neurnaseq, bg$NRN.LAM, range = 0, fun =lessthanz )
damid.tss.3(glirnaseq, bg$Glia.LAM, range = 0, fun =lessthanz )

damid.tss.3(neurnaseq, bg$NRN.PC, range = 0, fun =lessthanz )
damid.tss.3(glirnaseq, bg$Glia.PC, range = 0, fun =lessthanz )
