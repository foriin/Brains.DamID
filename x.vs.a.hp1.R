library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(reshape2)

rm(list = ls())

setwd(paste0("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
"DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/HP1.bedgraph.plus.Kc"))

hp1.bg <- dir(pattern = "HP1")
euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

# "Boundaries" of euchromatin

euc.gr <- GRanges(
  seqnames = Rle(euc),
  ranges = IRanges(c(1, 1600000, 1, 1, 1), end = c(22000000, 21146708, 22900000, 27905053, 22300000))
)

# no granges
bg.pl <- new.env()

for (i in hp1.bg){
  df <- fread(file.path(i), skip = 1, col.names = c("chr", "start", "end", "log2damid"))
  df$chr <- factor(df$chr)
  levels(df$chr) <- list(A = euc[1:4], X = "chrX")
  assign(sub("^[^.]+\\.([^.]+).*", "\\1", i), filter(df, !is.na(chr)), e = bg.pl)
}

boxplot(log2damid~ chr, data = bg.pl$BR, outline = F)

# domains

damid.fol <- paste0("~/IMG/Projects/",
"HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
"DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
"HP1.bedgraph.plus.Kc ")

domains <- new.env()
# hp1 domains
for (i in dir(damid.fol, pattern = "(HP1).*bed")){
  df <- fread(file.path(damid.fol, i), skip = 1, col.names = c("chr", "start", "end")) %>% dplyr::filter(chr %in%euc)
  assign(sub("\\.domains\\.bed", "", i),
         GenomicRanges::intersect(
           GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
           euc.gr),
         envir = domains)
}
rm(i)

# lam domains

dom.dir <- paste0("~/IMG/Projects/",
"HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
"DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
"BioHMM2.qn.full.PC.HMM3/")

for (i in dir(dom.dir, pattern = "(LAM).*bed")){
  df <- fread(file.path(dom.dir, i), skip = 1, col.names = c("chr", "start", "end")) %>% filter(chr %in%euc)
  assign(sub("\\.domains\\.bed", "", i),
         GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
         envir = domains)
}
rm(i)

# hp1 bedgraphs filtered by domains


bg <- new.env()

for (i in hp1.bg){
  df <- fread(i, skip = 1, col.names = c("chr", "start", "end", "log2damid")) %>% filter(chr %in%euc)
  assign(sub("^[^.]+\\.([^.]+).*", "\\1", i),
         subsetByOverlaps(GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end), damid = df$log2damid),
                          get(ls(envir = domains)[grepl(sub("^([^.]+)\\.([^.]+).*", "\\2.\\1", i), ls(e = domains))], e = domains),
                          minoverlap = 1001L),
         envir = bg)
}

in.lad.out.lad <- lapply(ls(envir = bg)[c(1,5,4,3,2)], function(nom){
 list(
   LADs = subsetByOverlaps(get(nom, envir = bg),
                    get(ls(envir = domains)[grepl(paste0(nom, ".LAM"), ls(envir = domains), ignore.case = T)], envir = domains),
                    minoverlap=1000L),
   intLADs = subsetByOverlaps(get(nom, envir = bg),
                    gaps(get(ls(envir = domains)[grepl(paste0(nom, ".LAM"), ls(envir = domains), ignore.case = T)], envir = domains)),
                    minoverlap=1000L)
 )
})

names(in.lad.out.lad) <- ls(envir = bg)[c(1,5,4,3,2)]

in.lad.df <- lapply(in.lad.out.lad, function(lst){
    chrom <- seqnames(lst[[1]])
    levels(chrom) <- list(A = euc[1:4], X = "chrX")
    data.frame(chr = chrom,
               start = start(lst[[1]]),
               end = end(lst[[1]]),
               damid = elementMetadata(lst[[1]])[,1]) %>% filter(!is.na(chr))
  })


out.lad.df <- lapply(in.lad.out.lad, function(lst){
  chrom <- seqnames(lst[[2]])
  levels(chrom) <- list(A = euc[1:4], X = "chrX")
  data.frame(chr = chrom,
             start = start(lst[[2]]),
             end = end(lst[[2]]),
             damid = elementMetadata(lst[[2]])[,1]) %>% filter(!is.na(chr))
})

bg.lst <- lapply(bg, function(gr){
  chrom <- seqnames(gr)
  levels(chrom) <- list(A = euc[1:4], X = "chrX")
  data.frame(chr = seqnames(gr),
             chrom = chrom,
             start = start(gr),
             end = end(gr),
             damid = elementMetadata(gr)[,1]) %>% 
    filter(!is.na(chrom))
})

bg.df <- Reduce(function(x, y) merge(x, y, by=c("chr", "start"), all.x = T, all.y = T), bg.lst)
levels(bg.df$chr) <- list(A = euc[1:4], X = "chrX")
bg.df <- bg.df[,c(1, 5, 8, 11, 14, 17)]
bg.df <- bg.df %>% setNames(c("chr", "BR", "NRN", "Kc167", "Glia", "FB"))
bg.melt <- melt(bg.df)

a <- melt(bg.df) %>% filter(variable == "damid") %>% droplevels()
pdf("a.vs.x.in.domains.pdf")
ggplot(a, aes(x = L1, y = value))+
  geom_boxplot(aes(fill = chr), position=position_dodge(1), outlier.shape = NA)+
  xlab("")+
  scale_y_continuous(limits = c(-3, 5))
dev.off()


b <- melt(in.lad.df) %>% filter(variable == "damid") %>% droplevels()

pdf("lads.pdf")
ggplot(b, aes(x = L1, y = value))+
  geom_boxplot(aes(fill = chr), position=position_dodge(1), outlier.shape = NA)+
  xlab("")+
  scale_y_continuous(limits = c(-3, 5.2))+
  ggtitle("LADs")
dev.off()

c <- melt(out.lad.df) %>% filter(variable == "damid") %>% droplevels()

pdf("interlads.pdf")
ggplot(c, aes(x = L1, y = value))+
  geom_boxplot(aes(fill = chr), position=position_dodge(1), outlier.shape = NA)+
  xlab("")+
  scale_y_continuous(limits = c(-2.5, 7))+
  ggtitle("interLADs")

dev.off()


# 
# boxplot(damid ~ chr, data = in.lad.out.lad.df$Kc167$intLADs, outline = F)
# 
# boxplot(damid ~ chr, data = bg.df$Kc167, outline = F)
# 
# sum(width(domains$BR.HP1))
# sum(width(bg$BR))
# sum(bg.pl$BR$end - bg.pl$BR$start)
# 
# bg.br.old <- fread("../Bedgraph/HP1.BR.m.sum.norm.qn.qn.bedgraph", skip = 1) %>% filter(V1 %in% euc)
# names(bg.br.old) <- c("chr", "start", "end", "sig")
# bg.br.old$chr <- factor(bg.br.old$chr)
# levels(bg.br.old$chr) <- list(A = euc[1:4], X = "chrX")
# bg.br.old <- filter(bg.br.old, !(is.na(chr)))
# 
# boxplot(sig ~ chr, data = bg.br.old, outline = F)


# P.values

sapply(bg.df, function(x) format(wilcox.test(damid ~ chr, x, alt = "l")$p.value, digits = 3))

sapply(in.lad.df, function(x) format(wilcox.test(damid ~ chr, x, alt = "l")$p.value, digits = 3))

sapply(out.lad.df, function(x) format(wilcox.test(damid ~ chr, x, alt = "l")$p.value, scientific = F, digits = 3))

pvs <- as.data.frame(t(sapply(list(bg.df,
                                 in.lad.df,
                                 out.lad.df), 
       function(y) {sapply(y, function(x) format(wilcox.test(damid ~ chr, x, alt = "l")$p.value, digits = 3))
}))) %>% setattr("row.names", c("Whole", "LADs", "interLADs"))

write.table(pvs, file = "x.vs.a.pvs.csv", sep = ";", quote = F, dec = ",", row.names = T, col.names = T)


# Ovaries

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/hp1_ov_som/")

hp1.som.ov <- fread("./23.10.17_Ovaries_HP1_1k/CSV/06.Dam.Normalized.DATA.csv") %>% filter(!(is.na(HP1.OV.wt.norm)))

hp1.som.ov$chr <- factor(paste0("chr", hp1.som.ov$chr))
levels(hp1.som.ov$chr) <- list(A = euc[1:4], X = "chrX")

hp1.som.ov <- hp1.som.ov %>% filter(!(is.na(chr)))

pdf("hp1.som.ov.pdf", width = 5)
  boxplot(HP1.OV.wt.norm ~ chr, data = hp1.som.ov, outline = F)
dev.off()
wilcox.test(HP1.OV.wt.norm ~ chr, data = hp1.som.ov, alt = "l")

