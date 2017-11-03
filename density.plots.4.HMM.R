library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)

rm(list = ls())


setwd(paste0("~/IMG/Projects/",
             "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
             "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
             "BioHMM2.qn.full.PC.HMM3/"))
load("densities.RData")


beds <- new.env()
for (i in dir()[grepl("bed", dir())]){
  assign(sub("\\.domains.bed", "", i), fread( i, skip = 1), envir = beds)
}
rm(i)

x <- sapply(names(beds), function(nym){
  df <- get(nym, envir = beds) %>% setNames(c("chr", "start", "end"))
  df$end - df$start
})

y <- melt(x) %>% setNames(c("width", "name"))


a <- ggplot(y[grepl("NRN", y$name),], aes(log10(width), col = name))+
  geom_density(adjust=1/2)+
  xlim(2, 6)+
  ggtitle("Neurons")+
  xlab(expression('log'[10]*'(domain size, bp)'))+
  scale_color_manual(labels = c("HP1", "LAM", "PC"), values = c("blue", "red", "orange"))

b <- ggplot(y[grepl("BR", y$name),], aes(log10(width), col = name))+
  geom_density()+
  xlim(2, 6)+
  ggtitle("Whole Brains")+
  xlab(expression('log'[10]*'(domain size, bp)'))+
  scale_color_manual(labels = c("HP1", "LAM", "PC"), values = c("blue", "red", "orange"))

c <- ggplot(y[grepl("Glia", y$name),], aes(log10(width), col = name))+
  geom_density()+
  xlim(2, 6)+
  ggtitle("Glia")+
  xlab(expression('log'[10]*'(domain size, bp)'))+
  scale_color_manual(labels = c("HP1", "LAM", "PC"), values = c("blue", "red", "orange"))

d <- ggplot(y[grepl("FB", y$name),], aes(log10(width), col = name))+
  geom_density()+
  xlim(2, 6)+
  ggtitle("Fatbody")+
  xlab(expression('log'[10]*'(domain size, bp)'))+
  scale_color_manual(labels = c("HP1", "LAM", "PC"), values = c("blue", "red", "orange"))



e <- ggplot(y[grepl("Kc", y$name),], aes(log10(width), col = name))+
  geom_density()+
  xlim(0, 6)+
  ggtitle("Kc167")+
  xlab(expression('log'[10]*'(domain size, bp)'))+
  scale_color_manual(labels = c("HP1", "LAM", "PC"), values = c("blue", "red", "orange"))

pdf("dens.plots.pdf", width = 12)
  grid.arrange(a, b, c, d, e, ncol = 3)
dev.off()

dom.stat <- y %>% group_by(name) %>% summarise(number = n(), median = median(width))

write.table(dom.stat, "domains.numbers.and.medians.csv", quote = F,
                        row.names = F, sep = ";", dec = ",")

a <- data.frame(prop.table(table(width = nrn$width, prot = nrn$name), 1))

ggplot(a, aes(x = width, y = Freq))+
  geom_line(aes(group = prot, color = prot))

ggplot(nrn, aes(x=log10(width))) + 
  geom_freqpoly(      # Histogram with density instead of count on y-axis
                 binwidth=0.4) +
  geom_density( alpha=.2)
