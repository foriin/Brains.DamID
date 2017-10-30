library(GenomicRanges)
library(data.table)
library(combinat)
library(eulerr)

rm(list = ls())

setwd(paste0("~/IMG/Projects/",
             "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/21.07.17_1kbins_NRN/",
             "sum_rep_qn/BioHMM2.qn.full.PC.HMM3/"))


beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)

prots <- unique(sub("[^.]+\\.([^.]+)", "\\1", names(beds)))
tiss <- unique(sub("^([^.]+)\\..*", "\\1", names(beds)))

comp_by_tiss <- lapply(tiss, function(tis){
  retl <- list()
  beds.t.i <- which(grepl(tis,  names(beds)))
  tis.grs <- lapply(beds.t.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    GRanges(seqnames = Rle(df$chr),
            ranges = IRanges(start = df$start, end = df$end))
  })
  names(tis.grs) <- sub("^.+\\.", "", names(beds)[beds.t.i])
  tis.grs <- tis.grs[order(names(tis.grs))]
 
  # Length of plain domains
  tis.grs.wid <- sapply(tis.grs, function(lst) sum(width(lst)))
 
  names(tis.grs.wid) <- names(tis.grs)
  retl[["dom.1x.width"]] <- tis.grs.wid
  
  # Width of all 2x intersects
  tis.grs.by2 <- combn(tis.grs, 2, simplify = F)
  
  tis.gr.width.by2 <- sapply(tis.grs.by2, function(lst){
    intrsct <- Reduce(GenomicRanges::intersect, lst)
    sum(width(reduce(intrsct)))
  })
  
  twb.names <- sapply(tis.grs.by2, function(comb){
    Reduce(function(x,y) paste(x, y, sep=".vs."), names(comb))
  })
  names(tis.gr.width.by2) <- twb.names
  retl[["int.dom.2x"]] <- (tis.gr.width.by2)
  
  # Width of 3x intersect
  tot.intrsct <- Reduce(GenomicRanges::intersect, tis.grs)
  retl[["int.dom.3x"]] <- sum(width(reduce(tot.intrsct)))
  
  # anti-domains
  tis.gaps <- lapply(tis.grs, gaps)
  gaps.intrsct <- Reduce(GenomicRanges::intersect, tis.gaps)
  retl[["antidomains.width"]] <- sum(width(reduce(gaps.intrsct)))
  return(retl)
})

names(comp_by_tiss) <- tiss

###################################

tiss.names <- data.frame(code = tiss, name =c("Neurons", "Whole brains", "Glia", "Kc167", "Fat bodies"))

x <- lapply(tiss, function(nym){
  lst <- comp_by_tiss[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  fit1 <- euler(c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x),
                  "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x),
                  "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x),
                  "HP1&Lam" = unname(intersect(hp1, lam) - tot.x),
                  "Pc&HP1" = unname(intersect(hp1, pc) - tot.x),
                  "Lam&Pc" = unname(intersect(pc, lam) - tot.x),
                  "Pc&Lam&HP1" = tot.x))
  fit1$original.values <- paste0(round(fit1$original.values/120381546*100, digits = 2), "%")
  # str(fit1)
  name <- tiss.names$name[tiss.names$code == nym]
  print(name)
  pl <- plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = name)
  pdf(paste0("~/", nym, ".pdf"))
    print(pl)
  dev.off
})


fit1 <- euler(c("HP1" = unname(comp_by_tiss[[1]][[1]][1] - comp_by_tiss[[1]][[2]][1] - comp_by_tiss[[1]][[2]][2] + comp_by_tiss[[1]][[3]]),
                "Lam" = unname(comp_by_tiss[[1]][[1]][2] - comp_by_tiss[[1]][[2]][3] - comp_by_tiss[[1]][[2]][1] + comp_by_tiss[[1]][[3]]),
                "Pc" = unname(comp_by_tiss[[1]][[1]][3] - comp_by_tiss[[1]][[2]][2] - comp_by_tiss[[1]][[2]][3] + comp_by_tiss[[1]][[3]]),
                "HP1&Lam" = unname(comp_by_tiss[[1]][[2]][1] - comp_by_tiss[[1]][[3]]),
                "Pc&HP1" = unname(comp_by_tiss[[1]][[2]][2] - comp_by_tiss[[1]][[3]]),
                "Lam&Pc" = unname(comp_by_tiss[[1]][[2]][3] - comp_by_tiss[[1]][[3]]),
                "Pc&Lam&HP1" = unname(comp_by_tiss[[1]][[3]])))
fit1$original.values <- paste0(round(fit1$original.values/120381546*100, digits = 2), "%")

for (i in 1:5){
  areas <- c("HP1" = unname(comp_by_tiss[[i]][[1]][1] - comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[2]][2] + comp_by_tiss[[i]][[3]]),
             "Lam" = unname(comp_by_tiss[[i]][[1]][2] - comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[2]][1] + comp_by_tiss[[i]][[3]]),
             "HP1&Lam" = unname(comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[3]]),
             "Pc" = unname(comp_by_tiss[[i]][[1]][3] - comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[2]][3] + comp_by_tiss[[i]][[3]]),
             "Pc&HP1" = unname(comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[3]]),
             "Lam&Pc" = unname(comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[3]]),
             "Pc&Lam&HP1" = unname(comp_by_tiss[[i]][[3]]))
  percents <- paste0(round(areas/120381546*100, digits = 2), "%")
  cat(tiss[i], "\n", names(areas), "\n", areas, "\n", percents, "\n", file = "venn.params.txt", append = T)
}




plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = "Neurons")


venn("100 + 110 + 101 + 111")
venn(5, ellipse = TRUE, zcolor="style", snames = c("BR", "Kc", "Glia", "FB", "Neur"))

x <- list(HP1 = c(0.2, 0.33, 0.05), LAM = c(0.1, 0.33, 0.05), Pc = c(0.2, 0.1, 0.05))
venn(x, zcolor = "style", ilabels = T)

set.seed(12345)
x <- as.data.frame(matrix(sample(0:1, 150, replace=TRUE), ncol=5))
venn(x)

venn.plot <- draw.quintuple.venn(
  area1 = 301,
  area2 = 321,
  area3 = 311,
  area4 = 321,
  area5 = 301,
  n12 = 188,
  n13 = 191,
  n14 = 184,
  n15 = 177,
  n23 = 194,
  n24 = 197,
  n25 = 190,
  n34 = 190,
  n35 = 173,
  n45 = 186,
  n123 = 112,
  n124 = 108,
  n125 = 108,
  n134 = 111,
  n135 = 104,
  n145 = 104,
  n234 = 111,
  n235 = 107,
  n245 = 110,
  n345 = 100,
  n1234 = 61,
  n1235 = 60,
  n1245 = 59,
  n1345 = 58,
  n2345 = 57,
  n12345 = 31,
  category = c("Neurons", "Glia", "Whole Brains", "Fat bodies", "Kc167"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
)

comp_by_tiss[[1]][1]/120381546*100
