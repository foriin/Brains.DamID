library("dplyr")
library("data.table")
library(GenomicRanges)
rm(list = ls())

sum_bed <- function(df){
  chr.len <- data.frame("chr" = c("2L", "2R", "3L", "3R", "4", "X", "2LHet",
                                  "2RHet", "3LHet", "3RHet", "XHet", "YHet"),
                        "len" = c(23011544, 21146708, 24543557, 27905053, 1351857,
                                  22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038))
  df %>% 
    mutate(len = end - start) %>% 
    group_by(chr) %>% 
    summarise(dom_len = sum(len)) %>% 
    rowwise() %>% 
    mutate(chrlen = chr.len$len[grepl(paste0("^", chr, "$"), chr.len$chr)])
}

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/02.04.17_all_samples_300/hmm3/")

beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.bed", "", i), fread(i), envir = beds)
}
rm(i)

chr.len <- data.frame("chr" = c("2L", "2R", "3L", "3R", "4", "X", "2LHet",
                                "2RHet", "3LHet", "3RHet", "XHet", "YHet"),
                      "len" = c(23011544, 21146708, 24543557, 27905053, 1351857,
                                22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038)
)

comparisons <- new.env()

for (i in 1:(length(ls(envir = beds)) - 1)){
  x <- get(ls(envir = beds)[i], envir = beds)
  names(x) <- c("chr", "start", "end")
  x.name <- ls(envir = beds)[i]
  lapply((i+1):length(ls(envir = beds)), function(k){
    compto <- get(ls(envir = beds)[k], envir = beds)
    names(compto) <- c("chr", "start", "end")
    compto.name <- ls(envir = beds)[k]
    gr1 <- GRanges(seqnames = Rle(x$chr),  ranges = IRanges(start = x$start, end = x$end))
    gr2 <- GRanges(seqnames = Rle(compto$chr),  ranges = IRanges(start = compto$start, end = compto$end))
    
    intersect.gr <- intersect(gr1, gr2)
    df_intersected <- data.frame(chr = seqnames(intersect.gr), start = start(intersect.gr), end = end(intersect.gr))
    ab <- merge(sum_bed(x)[,-3], sum_bed(compto)[,-3], by = "chr")
    y <- merge(ab, sum_bed(df_intersected), by = "chr") %>% 
      mutate("i/x" = round(dom_len / dom_len.x, digits=3),
             "i/y" = round(dom_len / dom_len.y, digits=3),
             coverage = round(dom_len / chrlen, digits=3))
    assign(paste0(x.name, ".vs.", compto.name), y, envir = comparisons)
    return()
      
  })
}

dir.create("Comparisons")
setwd("Comparisons/")

for (i in ls(envir = comparisons)){
  write.table(get(i, envir = comparisons), paste0(i, ".csv"), quote = F, sep = ";", dec = ",", row.names = F)
}


BR.HP1 %>% 
  setNames(c("chr", "start", "end")) %>% 
  mutate(len = end - start) %>% 
  group_by(chr) %>% 
  summarise(dom_len = sum(len)) %>% 
  rowwise() %>% 
  mutate(chrlen = chr.len$len[grepl(paste0("^", chr, "$"), chr.len$chr)]) %>% 
  mutate(ratio = dom_len / chrlen)

lapply(1:2, function(i){
  y <- get(ls(envir = beds)[i], envir = beds)
  name <- paste0(ls(envir = beds)[i], i)
  assign(name, y, envir = .GlobalEnv)
  return()
})
