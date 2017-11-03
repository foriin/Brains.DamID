# Generate venn diagrams for each tissue
# Euchromatin and heterochromatin are separated
# 4th chromosome is separated from rest

library(GenomicRanges)
library(data.table)
library(combinat)
library(eulerr)
library(dplyr)

rm(list = ls())

setwd(paste0("~/IMG/Projects/",
             "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
             "BioHMM2.qn.full.PC.HMM3/"))

# load workspace
load("venn3.RData")

euchr <- c("chr2L", "chr2R", "chr3L", "chr3R","chrX")
het <- c("chr2LHet", "chr2RHet", "chr3LHet", "chr3RHet","chrXHet")

# chrsize <- fread("~/IMG/data/dmel/Genome/dm3.genome")
# chrsize.h <- chrsize %>% filter(V1 %in% het)
# sum(chrsize.h$V2)


beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)

prots <- unique(sub("[^.]+\\.([^.]+)", "\\1", names(beds)))
tiss <- unique(sub("^([^.]+)\\..*", "\\1", names(beds)))
tissnames <- c("Neurons", "Whole Brains", "Glia", "Kc167", "Fat Bodies")

# Euchromatin

comp_by_tiss <- lapply(tiss, function(tis){
  retl <- list()
  beds.t.i <- which(grepl(tis,  names(beds)))
  tis.grs <- lapply(beds.t.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr %in% euchr)
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

tiss.names <- data.frame(code = tiss, name = tissnames)

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
  fit1$original.values <- paste0(round(fit1$original.values/119029689*100, digits = 2), "%")
  # str(fit1)
  name <- tiss.names$name[tiss.names$code == nym]
  print(name)
  pl <- plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = name)
  pdf(paste0("~/", nym, ".pdf"))
    print(pl)
  dev.off
})

info.list <- sapply(tiss, function(nym){
  lst <- comp_by_tiss[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  diviser = c(hp1[1], lam[1], pc[1])
  round(sapply(diviser, function(div){
    c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x)/div,
      "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x)/div,
      "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x)/div,
      "HP1&Lam" = unname(intersect(hp1, lam) - tot.x)/div,
      "Pc&HP1" = unname(intersect(hp1, pc) - tot.x)/div,
      "Lam&Pc" = unname(intersect(pc, lam) - tot.x)/div,
      "Pc&Lam&HP1" = tot.x/div)
  })*100, digits = 2)
  
}, simplify = F)

info.list <- lapply(names(info.list), function(nym){
  df <- as.data.frame(info.list[[nym]])
  colnames(df) <- paste0(nym, ".", colnames(df))
  df
})

info.table <- Reduce("cbind", info.list)
type_of_i <- row.names(info.table)

info.table <- cbind(type_of_i, sapply(names(info.table), function(nom){
  y <- info.table[[nom]]
  y[!grepl(sub("^[^.]+\\.(.*)$", "\\1", nom), row.names(info.table), ignore.case = T)] <- NA
  y
}))

# apply(info.table, 2, sum, na.rm = T)

write.table(info.table, "euc.intersections.info.csv", quote = F, sep = ";", dec = ",",
            row.names = F)


fit1 <- euler(c("HP1" = unname(comp_by_tiss[[1]][[1]][1] - comp_by_tiss[[1]][[2]][1] - comp_by_tiss[[1]][[2]][2] + comp_by_tiss[[1]][[3]]),
                "Lam" = unname(comp_by_tiss[[1]][[1]][2] - comp_by_tiss[[1]][[2]][3] - comp_by_tiss[[1]][[2]][1] + comp_by_tiss[[1]][[3]]),
                "Pc" = unname(comp_by_tiss[[1]][[1]][3] - comp_by_tiss[[1]][[2]][2] - comp_by_tiss[[1]][[2]][3] + comp_by_tiss[[1]][[3]]),
                "HP1&Lam" = unname(comp_by_tiss[[1]][[2]][1] - comp_by_tiss[[1]][[3]]),
                "Pc&HP1" = unname(comp_by_tiss[[1]][[2]][2] - comp_by_tiss[[1]][[3]]),
                "Lam&Pc" = unname(comp_by_tiss[[1]][[2]][3] - comp_by_tiss[[1]][[3]]),
                "Pc&Lam&HP1" = unname(comp_by_tiss[[1]][[3]])))
fit1$original.values <- paste0(round(fit1$original.values/119029689*100, digits = 2), "%")

for (i in 1:5){
  fit <- euler(c("HP1" = unname(comp_by_tiss[[i]][[1]][1] - comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[2]][2] + comp_by_tiss[[i]][[3]]),
             "Lam" = unname(comp_by_tiss[[i]][[1]][2] - comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[2]][1] + comp_by_tiss[[i]][[3]]),
             "HP1&Lam" = unname(comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[3]]),
             "Pc" = unname(comp_by_tiss[[i]][[1]][3] - comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[2]][3] + comp_by_tiss[[i]][[3]]),
             "Pc&HP1" = unname(comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[3]]),
             "Lam&Pc" = unname(comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[3]]),
             "Pc&Lam&HP1" = unname(comp_by_tiss[[i]][[3]])))
  fit$original.values <- paste0(round(fit$original.values/119029689*100, digits = 2), "%")
  pdf(file = paste0(names(comp_by_tiss)[i], ".venn.pdf"))
    plot(fit, fill = c("blue", "red", "khaki"), counts = T, main = tissnames[i])
  dev.off()
}

for (i in 1:5){
  areas <- c("HP1" = unname(comp_by_tiss[[i]][[1]][1] - comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[2]][2] + comp_by_tiss[[i]][[3]]),
             "Lam" = unname(comp_by_tiss[[i]][[1]][2] - comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[2]][1] + comp_by_tiss[[i]][[3]]),
             "HP1&Lam" = unname(comp_by_tiss[[i]][[2]][1] - comp_by_tiss[[i]][[3]]),
             "Pc" = unname(comp_by_tiss[[i]][[1]][3] - comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[2]][3] + comp_by_tiss[[i]][[3]]),
             "Pc&HP1" = unname(comp_by_tiss[[i]][[2]][2] - comp_by_tiss[[i]][[3]]),
             "Lam&Pc" = unname(comp_by_tiss[[i]][[2]][3] - comp_by_tiss[[i]][[3]]),
             "Pc&Lam&HP1" = unname(comp_by_tiss[[i]][[3]]))
  percents <- paste0(round(areas/119029689*100, digits = 2), "%")
  cat(tiss[i], "\n", names(areas), "\n", areas, "\n", percents, "\n", file = "venn.params.txt", append = T)
}




plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = "Neurons")

# Only for X chromosome

comp_by_tiss.X <- lapply(tiss, function(tis){
  retl <- list()
  beds.t.i <- which(grepl(tis,  names(beds)))
  tis.grs <- lapply(beds.t.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr == "chrX")
    GRanges(seqnames = Rle(df$chr),
            ranges = IRanges(start = df$start, end = df$end))
  })
  print(tis)
  print(sapply(tis.grs, function(x) sum(width(x))))
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

names(comp_by_tiss.X) <- tiss
tiss.names <- data.frame(code = tiss, name = tissnames)

x <- lapply(tiss, function(nym){
  lst <- comp_by_tiss.X[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  test <- c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x),
            "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x),
            "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x),
            "HP1&Lam" = unname(intersect(hp1, lam) - tot.x),
            "Pc&HP1" = unname(intersect(hp1, pc) - tot.x),
            "Lam&Pc" = unname(intersect(pc, lam) - tot.x),
            "Pc&Lam&HP1" = tot.x)
  # print(c(hp1, intersect(hp1, lam), tot.x))
  # print(test)
  fit1 <- euler(test)
  fit1$original.values <- paste0(round(fit1$original.values/22422827*100, digits = 2), "%")
  # str(fit1)
  name <- tiss.names$name[tiss.names$code == nym]
  print(name)
  fit1
  # pl <- plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = name)
  # pdf(paste0("~/", nym, ".pdf"))
  # print(pl)
  # dev.off
})

pdf("nrn.X.pdf")
plot(x[[1]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("br.X.pdf")
plot(x[[2]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("glia.X.pdf")
plot(x[[3]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("kc.X.pdf")
plot(x[[4]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("fb.X.pdf")
plot(x[[5]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

info.list <- sapply(tiss, function(nym){
  lst <- comp_by_tiss.X[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  diviser = c(hp1[1], lam[1], pc[1])
  round(sapply(diviser, function(div){
    c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x)/div,
    "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x)/div,
    "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x)/div,
    "HP1&Lam" = unname(intersect(hp1, lam) - tot.x)/div,
    "Pc&HP1" = unname(intersect(hp1, pc) - tot.x)/div,
    "Lam&Pc" = unname(intersect(pc, lam) - tot.x)/div,
    "Pc&Lam&HP1" = tot.x/div)
  })*100, digits = 2)
  
}, simplify = F)

info.list <- lapply(names(info.list), function(nym){
  df <- as.data.frame(info.list[[nym]])
  colnames(df) <- paste0(nym, ".", colnames(df))
  df
})

info.table <- Reduce("cbind", info.list)
type_of_i <- row.names(info.table)

info.table <- cbind(type_of_i, sapply(names(info.table), function(nom){
  y <- info.table[[nom]]
  y[!grepl(sub("^[^.]+\\.(.*)$", "\\1", nom), row.names(info.table), ignore.case = T)] <- NA
  y
}))

# apply(info.table, 2, sum, na.rm = T)

write.table(info.table, "X.intersections.info.csv", quote = F, sep = ";", dec = ",",
            row.names = F)


# Heterochromatin

tiss <- tiss[c(1:3, 5)]

comp_by_tiss.h <- lapply(tiss, function(tis){
  retl <- list()
  beds.t.i <- which(grepl(tis,  names(beds)))
  tis.grs <- lapply(beds.t.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr %in% het)
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

names(comp_by_tiss.h) <- tiss

tiss.names <- data.frame(code = tiss, name = tissnames[c(1:3, 5)])

x <- lapply(tiss[c(1:3, 5)], function(nym){
  lst <- comp_by_tiss.h[[nym]]
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
  fit1$original.values <- paste0(round(fit1$original.values/8934743*100, digits = 2), "%")
  # str(fit1)
  name <- tiss.names$name[tiss.names$code == nym]
  print(name)
  fit1
  # pl <- plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = name)
  # pdf(paste0("~/", nym, ".pdf"))
  # print(pl)
  # dev.off
})

pdf("nrn.het.pdf")
plot(x[[1]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("br.het.pdf")
plot(x[[2]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("glia.het.pdf")
plot(x[[3]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("fb.het.pdf")
plot(x[[4]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

sapply(tiss, function(nym){
  lst <- comp_by_tiss.h[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
 c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x),
                  "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x),
                  "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x),
                  "HP1&Lam" = unname(intersect(hp1, lam) - tot.x),
                  "Pc&HP1" = unname(intersect(hp1, pc) - tot.x),
                  "Lam&Pc" = unname(intersect(pc, lam) - tot.x),
                  "Pc&Lam&HP1" = tot.x)

})

info.list <- sapply(tiss, function(nym){
  lst <- comp_by_tiss.h[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  diviser = c(hp1[1], lam[1], pc[1])
  round(sapply(diviser, function(div){
    c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x)/div,
      "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x)/div,
      "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x)/div,
      "HP1&Lam" = unname(intersect(hp1, lam) - tot.x)/div,
      "Pc&HP1" = unname(intersect(hp1, pc) - tot.x)/div,
      "Lam&Pc" = unname(intersect(pc, lam) - tot.x)/div,
      "Pc&Lam&HP1" = tot.x/div)
  })*100, digits = 2)
  
}, simplify = F)

info.list <- lapply(names(info.list), function(nym){
  df <- as.data.frame(info.list[[nym]])
  colnames(df) <- paste0(nym, ".", colnames(df))
  df
})

info.table <- Reduce("cbind", info.list)
type_of_i <- row.names(info.table)

info.table <- cbind(type_of_i, sapply(names(info.table), function(nom){
  y <- info.table[[nom]]
  y[!grepl(sub("^[^.]+\\.(.*)$", "\\1", nom), row.names(info.table), ignore.case = T)] <- NA
  y
}))

# apply(info.table, 2, sum, na.rm = T)

write.table(info.table, "Het.intersections.info.csv", quote = F, sep = ";", dec = ",",
            row.names = F)


# 4th chromosome

setwd(paste0("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
"DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
"BioHMM2.qn.full.PC.HMM3/4th_chrom/BioHMM.qn/"))

tiss <- unique(sub("^([^.]+)\\..*", "\\1", names(beds)))
tissnames <- c("Neurons", "Whole Brains", "Glia", "Kc167", "Fat Bodies")

beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.bed", "", i), fread(i, skip = 1), envir = beds)
}
rm(i)


comp_by_tiss.4 <- lapply(tiss, function(tis){
  retl <- list()
  beds.t.i <- which(grepl(tis,  names(beds)))
  tis.grs <- lapply(beds.t.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr == "chr4" | chr == "4")
    GRanges(seqnames = Rle(df$chr),
            ranges = IRanges(start = df$start, end = df$end))
  })
  print(tis)
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

names(comp_by_tiss.4) <- tiss
tiss.names <- data.frame(code = tiss, name = tissnames)

x <- lapply(tiss, function(nym){
  lst <- comp_by_tiss.4[[nym]]
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
  fit1$original.values <- paste0(round(fit1$original.values/1351857*100, digits = 2), "%")
  # str(fit1)
  name <- tiss.names$name[tiss.names$code == nym]
  print(name)
  fit1
  # pl <- plot(fit1, fill = c("blue", "red", "khaki"), counts = T, main = name)
  # pdf(paste0("~/", nym, ".pdf"))
  # print(pl)
  # dev.off
})

pdf("nrn.4.pdf")
plot(x[[1]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("br.4.pdf")
plot(x[[2]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("glia.4.pdf")
plot(x[[3]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("kc.4.pdf")
plot(x[[4]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

pdf("fb.4.pdf")
plot(x[[5]], fill = c("blue", "red", "khaki"), counts = T)
dev.off()

info.list <- sapply(tiss, function(nym){
  lst <- comp_by_tiss.4[[nym]]
  hp1 <- Reduce(c, sapply(lst, function(x) x[grepl("HP1", names(x))]))
  lam <- Reduce(c, sapply(lst, function(x) x[grepl("LAM", names(x))]))
  pc <- Reduce(c, sapply(lst, function(x) x[grepl("PC", names(x))]))
  tot.x <- lst[[3]]
  diviser = c(hp1[1], lam[1], pc[1])
  round(sapply(diviser, function(div){
    c("HP1" = unname(hp1[1] - hp1[2] - hp1[3] + tot.x)/div,
      "Lam" = unname(lam[1] - lam[2] - lam[3] + tot.x)/div,
      "Pc" = unname(pc[1] - pc[2] - pc[3] + tot.x)/div,
      "HP1&Lam" = unname(intersect(hp1, lam) - tot.x)/div,
      "Pc&HP1" = unname(intersect(hp1, pc) - tot.x)/div,
      "Lam&Pc" = unname(intersect(pc, lam) - tot.x)/div,
      "Pc&Lam&HP1" = tot.x/div)
  })*100, digits = 2)
  
}, simplify = F)

info.list <- lapply(names(info.list), function(nym){
  df <- as.data.frame(info.list[[nym]])
  colnames(df) <- paste0(nym, ".", colnames(df))
  df
})

info.table <- Reduce("cbind", info.list)
type_of_i <- row.names(info.table)

info.table <- cbind(type_of_i, as.data.frame(sapply(names(info.table), function(nom){
  y <- info.table[[nom]]
  y[!grepl(sub("^[^.]+\\.(.*)$", "\\1", nom), row.names(info.table), ignore.case = T)] <- NA
  y
})))

# apply(info.table, 2, sum, na.rm = T)

write.table(info.table, "chr4.intersections.info.csv", quote = F, sep = ";", dec = ",",
            row.names = F)


