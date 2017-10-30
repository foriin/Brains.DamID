rm(list=ls())
# install.packages("exactRankTests")
library(exactRankTests)
library(data.table)
library(reshape2)
library(ggplot2)

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/HP1_on_X/06.09.17.NRN_GLIA_FATB_HP1_Lam_Pc_plus_Kc167_HP1_1k/Bedgraph/")

euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")


utest.box <- function(dataname){
  data <- fread(dataname, skip = 1) %>% setNames(c("chr", "start", "end",
                                                   sub("^[^.]+\\.([^.]+)\\..*$", "\\1", dataname))) %>%
    filter(chr %in% euc) %>% mutate(chr = as.factor(chr), chrom = chr) %>% select(-end)
  levels(data$chrom) <- list(A = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4"), X="chrX")
  
  data
  # cat(sub("(^[^.]+\\.[^.]+)\\..*$", "\\1", dataname), "\t", wilcox.test(aut, x, alt = "l", correct = T)$p.value, "\n",
  #     file = "A_vs_X_wilcox.plusKc.txt", append = T)
  # pdf(file.path("./PDF", sub("(^[^.]+\\.[^.]+)\\..*$", "\\1.boxplot.pdf", dataname)))
  #   boxplot(aut, x, outline = F, names = c("A", "X"), ylim=c(-6, 6))
  # dev.off()
}

chr.dat <- lapply(dir(pattern = "bedgraph"), utest.box)
chr.dat <- chr.dat[c(1,5,3,2,4)]

all.tis <- Reduce(function(x,y) merge(x, y, by = c("chr", "chrom", "start"), all = T), chr.dat) %>% 
  select(-chr, -start)

sapply(all.tis[2:6], function(col){
  wilcox.test(col[all.tis[,1] == "A"], col[all.tis[,1] == "X"], alt = "l")$p.value
})

all.tis.m <- melt(all.tis, id.var = "chrom", na.rm = T)

p <- ggplot(all.tis.m, aes(x=variable, y=value)) + 
  theme_bw()+
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))+
  geom_boxplot(outlier.shape = NA, aes(fill=chrom), position=position_dodge(width=0.8)) +
  ggtitle("HP1 on X and autosomes")

pdf("HP1_on_X_vs_A.bxp.pdf")
  p
dev.off()

