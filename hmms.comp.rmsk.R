library(data.table)
library(dplyr)
library(GenomicRanges)

setwd("~/IMG/data/dmel/repeats/")

ltrs <- fread("lines.ltrs.merged.bed")

names(ltrs) <- c("chr", "start", "end")

len.tran <- ltrs %>% group_by(chr) %>% 
  mutate(len = end - start) %>% 
  summarize(len = sum(len))

chr.len <- data.frame("chr" = c("2L", "2R", "3L", "3R", "4", "X", "2LHet",
                                "2RHet", "3LHet", "3RHet", "XHet", "YHet"),
                      "len" = c(23011544, 21146708, 24543557, 27905053, 1351857,
                                22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038))
chr.len$chr <- paste0("chr", chr.len$chr)

chrlen <- merge(chr.len, len.tran, by = "chr") %>% 
  transmute(chr = chr, len = len.x - len.y)

sum_bed <- function(df, chr.len = chrlen){
  
  df %>% 
    mutate(len = end - start) %>% 
    group_by(chr) %>% 
    summarise(dom_len = sum(len)) %>% 
    rowwise() %>% 
    mutate(chrlen = chr.len$len[grepl(paste0("^", chr, "$"), chr.len$chr)])
}


setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/23.06.17_with_quantnorm/BioHMM.ps/rmsk/")

beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.pseudo\\.domains\\.rmsk\\.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)



comparisons <- new.env()

for (i in 1:(length(ls(envir = beds)) - 1)){
  x <- get(ls(envir = beds)[i], envir = beds)
  names(x) <- c("chr", "start", "end")
  x.name <- ls(envir = beds)[i]
  lapply((i+1):length(ls(envir = beds)), function(k){
    compto <- get(ls(envir = beds)[k], envir = beds)
    names(compto) <- c("chr", "start", "end")
    compto.name <- ls(envir = beds)[k]
    print(compto.name)
    gr1 <- GRanges(seqnames = Rle(x$chr),  ranges = IRanges(start = x$start, end = x$end))
    gr2 <- GRanges(seqnames = Rle(compto$chr),  ranges = IRanges(start = compto$start, end = compto$end))
    
    intersect.gr <- intersect(gr1, gr2)
    df_intersected <- data.frame(chr = seqnames(intersect.gr), start = start(intersect.gr), end = end(intersect.gr))
    ab <- merge(sum_bed(x), sum_bed(compto), by = c("chr", "chrlen"), all.x = T, all.y = T)
    y <- merge(ab, sum_bed(df_intersected), by = c("chr", "chrlen"), all.x = T, all.y = T) %>% 
      mutate("i/x" = round(dom_len / dom_len.x, digits=3),
             "i/y" = round(dom_len / dom_len.y, digits=3),
             coverage.x = round(dom_len.x / chrlen, digits=3),
             coverage.y = round(dom_len.y / chrlen, digits=3),
             coverage.i = round(dom_len / chrlen, digits=3)) %>% 
      dplyr::rename(dom_len.i = dom_len) #bitch-ass s4vectors package hid dplyr rename function
    assign(paste0(x.name, ".vs.", compto.name), y, envir = comparisons)
    return()
    
  })
}

dir.create("Comparisons")
setwd("Comparisons/")

for (i in ls(envir = comparisons)){
  
  sumry.eu = get(i, envir = comparisons) %>% 
    filter(!grepl("Het|4", chr)) %>% 
    summarise_each(funs(sum), dom_len.x, dom_len.y, dom_len.i, chrlen) %>% 
    mutate("i/x total" = round(dom_len.i / dom_len.x, digits = 3), "i/y total" = round(dom_len.i / dom_len.y, digits = 3),
           "coverage x" = round(dom_len.x / chrlen, digits = 3),
           "coverage y" = round(dom_len.y / chrlen, digits = 3),
           "coverage i" = round(dom_len.i / chrlen, digits = 3))
  sumry.het = get(i, envir = comparisons) %>% 
    filter(grepl("Het|4", chr)) %>% 
    summarise_each(funs(sum(., na.rm = TRUE)), dom_len.x, dom_len.y, dom_len.i, chrlen) %>% 
    mutate("i/x total" = round(dom_len.i / dom_len.x, digits = 3), "i/y total" = round(dom_len.i / dom_len.y, digits = 3),
           "coverage x" = round(dom_len.x / chrlen, digits = 3),
           "coverage y" = round(dom_len.y / chrlen, digits = 3),
           "coverage i" = round(dom_len.i / chrlen, digits = 3))
  write.table(get(i, envir = comparisons), paste0(i, ".csv"), quote = F, sep = ";", dec = ",", row.names = F)
  cat("Euchromatin summary\n", file = paste0(i, ".csv"), append = T)
  write.table(sumry.eu, paste0(i, ".csv"), quote = F, sep = ";", dec = ",", row.names = F, append = T)
  cat("Heterochromatin summary\n", file = paste0(i, ".csv"), append = T)
  write.table(sumry.het, paste0(i, ".csv"), quote = F, sep = ";", dec = ",", row.names = F, append = T)
}

