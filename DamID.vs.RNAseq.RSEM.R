library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(ggsignif)


rm(list=ls())

euc <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/")

# Genes annotation
genes <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf") %>% 
  select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('.*gene_id "(FBgn[0-9]+)";.*', '\\1', attr), gene_name = sub('.*gene_name "([^;]+)";.*', '\\1', attr),
         tss = ifelse(strand == "+", start, end), chr = paste0("chr", chr)) %>% 
  select(-attr)%>% filter(chr %in% euc)

# Load rsem results for each tissue
br.rsem <- fread("BR_rsem.genes.results") %>% mutate(id = gsub("^([^_]+)_.*", "\\1", gene_id),
                                                     gene_name = gsub("^[^_]+_(.*)$", "\\1", gene_id)) %>% select(-gene_id) %>% 
  merge(., genes, by = c("id", "gene_name"))
fb.rsem <- fread("FB_rsem.genes.results") %>% mutate(id = gsub("^([^_]+)_.*", "\\1", gene_id),
                                                    gene_name = gsub("^[^_]+_(.*)$", "\\1", gene_id)) %>% select(-gene_id) %>% 
  merge(., genes, by = c("id", "gene_name"))

kc.rsem <- fread("KC_rsem.genes.results") %>% mutate(id = gsub("^([^_]+)_.*", "\\1", gene_id),
                                                       gene_name = gsub("^[^_]+_(.*)$", "\\1", gene_id)) %>% select(-gene_id) %>% 
  merge(., genes, by = c("id", "gene_name"))

# write.table(br.rsem, "Brains_expression_rsem.csv", quote = F,
#             row.names = F, sep = ";", dec = ",")
# 
# write.table(fb.rsem, "FatBodies_expression_rsem.csv", quote = F,
#             row.names = F, sep = ";", dec = ",")
# write.table(kc.rsem, "Kc167_expression_rsem.csv", quote = F,
#             row.names = F, sep = ";", dec = ",")

# Generate IRanges object for tss or whole genes
tss.gr <- GRanges(seqnames = Rle(genes$chr),  ranges = IRanges(start = genes$tss, width = 1, names = genes$id))
genes.gr <- GRanges(seqnames = Rle(genes$chr), ranges = IRanges(start = genes$start, end = genes$end, names = genes$id))



###################################
# DamID data
damid.fol <- "~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/BioHMM2.qn.full.PC.HMM3/"
domains <- new.env()
for (i in dir(damid.fol)[grepl("BR|FB|Kc167", dir(damid.fol))]){
   df <- fread(file.path(damid.fol, i), skip = 1, col.names = c("chr", "start", "end")) %>% filter(chr %in%euc)
   assign(sub("\\.domains\\.bed", "", i),
          GRanges(seqnames = Rle(df$chr),  ranges = IRanges(start = df$start, end = df$end)),
          envir = domains)
}
rm(i)

beds <- new.env()
for (i in dir(damid.fol)[grepl("bed", dir(damid.fol))]){
  assign(sub("\\.domains.bed", "", i), fread(file.path(damid.fol, i), skip = 1), envir = beds)
}
rm(i)

# I. Expression of genes located in both Lam and HP1 or Lam and Pc
#   or Lam - HP1 - Pc domains in brains and fatbodies

domains$BR.HP1.x.Lam <- GenomicRanges::intersect(domains$BR.HP1, domains$BR.LAM)
domains$BR.Pc.x.Lam <- GenomicRanges::intersect(domains$BR.PC, domains$BR.LAM)
domains$BR.SoloLamino <- GenomicRanges::setdiff(domains$BR.LAM, domains$BR.PC)
domains$BR.SoloLamino <- GenomicRanges::setdiff(domains$BR.SoloLamino, domains$BR.HP1)

domains$Kc167.SoloHP1 <- GenomicRanges::setdiff(domains$Kc167.HP1, domains$Kc167.LAM)
domains$Kc167.SoloHP1 <- GenomicRanges::setdiff(domains$Kc167.SoloHP1, domains$Kc167.PC)

br.damid.exp <- cbind(br.rsem, sapply(ls(pattern = "BR", envir = domains), function(grname){
  gr <- get(grname, envir = domains)
  in_domain <- subsetByOverlaps(tss.gr, gr)@ranges@NAMES
  ifelse(br.rsem$id %in% in_domain, 1, 0)
}))

# write.table(br.damid.exp, "br.rsem.domains.csv", quote = F, dec = ",",
#             sep = ";", row.names = F)



boxplot(br.damid.exp$TPM[br.damid.exp$BR.PC == 1 & br.damid.exp$TPM > 0], br.damid.exp$TPM[br.damid.exp$BR.PC == 0 & br.damid.exp$TPM > 0], outline = F)
boxplot(br.damid.exp$TPM[br.damid.exp$BR.LAM == 1 & br.damid.exp$TPM > 0], br.damid.exp$TPM[br.damid.exp$BR.LAM == 0 & br.damid.exp$TPM > 0], outline = F)
boxplot(br.damid.exp$TPM[br.damid.exp$BR.HP1 == 1 & br.damid.exp$TPM > 0], br.damid.exp$TPM[br.damid.exp$BR.HP1 == 0 & br.damid.exp$TPM > 0], outline = F)

# br.4.plot <- rbind(data.frame(tpm = br.damid.exp$TPM[br.damid.exp$BR.HP1.x.Lam == 1 ], type = "HP1.x.LAM"),
#                               data.frame(tpm = br.damid.exp$TPM[br.damid.exp$BR.SoloLamino == 1 ], type = "sLAM"),
#                               data.frame(tpm = br.damid.exp$TPM[br.damid.exp$BR.Pc.x.Lam == 1], type = "PC.x.LAM"))
# write.table(br.4.plot, "Brains.tpm.domains.for.mann.csv", quote = F,
#             row.names = F, sep = ";", dec = ",")

pdf("brains.lam.comb.pdf", width = 5.5)
  boxplot(br.damid.exp$TPM[br.damid.exp$BR.HP1.x.Lam == 1 & br.damid.exp$TPM > 0],
          br.damid.exp$TPM[br.damid.exp$BR.SoloLamino == 1 & br.damid.exp$TPM > 0],
          br.damid.exp$TPM[br.damid.exp$BR.Pc.x.Lam == 1 & br.damid.exp$TPM > 0],
          names = c("HP1 x Lam", "Lam - (HP1 + Pc)", "Pc x Lam"),
          main = "Brains", outline = F)
  text(1.3, 2.7, length(br.damid.exp$TPM[br.damid.exp$BR.HP1.x.Lam == 1 & br.damid.exp$TPM > 0]))
  text(2.3, 6, length(br.damid.exp$TPM[br.damid.exp$BR.SoloLamino == 1 & br.damid.exp$TPM > 0]))
  text(3.3, 2.9, length(br.damid.exp$TPM[br.damid.exp$BR.Pc.x.Lam == 1 & br.damid.exp$TPM > 0]))
dev.off()

wilcox.test(br.damid.exp$TPM[br.damid.exp$BR.HP1.x.Lam == 1 & br.damid.exp$TPM > 0],
            br.damid.exp$TPM[br.damid.exp$BR.SoloLamino == 1 & br.damid.exp$TPM > 0],
            alt = "l")

wilcox.test(br.damid.exp$TPM[br.damid.exp$BR.Pc.x.Lam == 1 & br.damid.exp$TPM > 0 & br.damid.exp$TPM > 0],
            br.damid.exp$TPM[br.damid.exp$BR.SoloLamino == 1 & br.damid.exp$TPM > 0 & br.damid.exp$TPM > 0],
            alt = "l")


domains$FB.HP1.x.Lam <- GenomicRanges::intersect(domains$FB.HP1, domains$FB.LAM)
domains$FB.Pc.x.Lam <- GenomicRanges::intersect(domains$FB.PC, domains$FB.LAM)
domains$FB.SoloLamino <- GenomicRanges::setdiff(domains$FB.LAM, domains$FB.PC)
domains$FB.SoloLamino <- GenomicRanges::setdiff(domains$FB.SoloLamino, domains$FB.HP1)


fb.damid.exp <- cbind(fb.rsem, sapply(ls(pattern = "FB", envir = domains), function(grname){
  gr <- get(grname, envir = domains)
  in_domain <- subsetByOverlaps(tss.gr, gr)@ranges@NAMES
  ifelse(fb.rsem$id %in% in_domain, 1, 0)
}))

# write.table(fb.damid.exp, "fb.rsem.domains.csv", quote = F, dec = ",",
#             sep = ";", row.names = F)

boxplot(fb.damid.exp$TPM[fb.damid.exp$FB.PC == 1 & fb.damid.exp$TPM > 0], fb.damid.exp$TPM[fb.damid.exp$FB.PC == 0 & fb.damid.exp$TPM > 0], outline = F)
boxplot(fb.damid.exp$TPM[fb.damid.exp$FB.LAM == 1 & fb.damid.exp$TPM > 0], fb.damid.exp$TPM[fb.damid.exp$FB.LAM == 0 & fb.damid.exp$TPM > 0], outline = F)
boxplot(fb.damid.exp$TPM[fb.damid.exp$FB.HP1 == 1 & fb.damid.exp$TPM > 0], fb.damid.exp$TPM[fb.damid.exp$FB.HP1 == 0 & fb.damid.exp$TPM > 0], outline = F)

# fb.4.plot <- rbind(data.frame(tpm = fb.damid.exp$TPM[fb.damid.exp$FB.HP1.x.Lam == 1 ], type = "HP1.x.LAM"),
#                    data.frame(tpm = fb.damid.exp$TPM[fb.damid.exp$FB.SoloLamino == 1 ], type = "sLAM"),
#                    data.frame(tpm = fb.damid.exp$TPM[fb.damid.exp$FB.Pc.x.Lam == 1], type = "PC.x.LAM"))
# write.table(fb.4.plot, "Fatbodies.tpm.domains.for.mann.csv", quote = F,
#             row.names = F, sep = ";", dec = ",")

pdf("fatbodies.lam.comb.pdf", width = 5.5)
  boxplot(fb.damid.exp$TPM[fb.damid.exp$FB.HP1.x.Lam == 1 ],
          fb.damid.exp$TPM[fb.damid.exp$FB.SoloLamino == 1 ],
          names = c("HP1 x Lam", "Lam - (HP1 + Pc)", "Pc x Lam"),
          main = "Fat Bodies", outline = F)
  text(1.3, 0.17, length(fb.damid.exp$TPM[fb.damid.exp$FB.HP1.x.Lam == 1 & fb.damid.exp$TPM > 0]))
  text(2.3, 0.21, length(fb.damid.exp$TPM[fb.damid.exp$FB.SoloLamino == 1 & fb.damid.exp$TPM > 0]))
  text(3.3, 0.24, length(fb.damid.exp$TPM[fb.damid.exp$FB.Pc.x.Lam == 1 & fb.damid.exp$TPM > 0]))
dev.off()

wilcox.test(fb.damid.exp$TPM[fb.damid.exp$FB.HP1.x.Lam == 1 ],
            fb.damid.exp$TPM[fb.damid.exp$FB.SoloLamino == 1 ],
            alt = "l")

wilcox.test(fb.damid.exp$TPM[fb.damid.exp$FB.Pc.x.Lam == 1  ],
            fb.damid.exp$TPM[fb.damid.exp$FB.SoloLamino == 1  ],
            alt = "g")

# II. Expression of genes in conservative across all tissues HP1 domains

prots <- unique(sub("[^.]+\\.([^.]+)$", "\\1", names(beds)))

totals.and.leftouts <- lapply(prots, function(prot){
  retl <- list()
  beds.p.i <- which(grepl(prot,  names(beds)))
  prots.grs <- lapply(beds.p.i, function(x){
    df <- get(names(beds)[x], envir = beds)
    names(df) <- c("chr", "start", "end")
    df <- df %>% filter(chr %in% euc)
    GRanges(seqnames = Rle(df$chr),
            ranges = IRanges(start = df$start, end = df$end))
  })
  names(prots.grs) <- sub("\\..+$", "", names(beds)[beds.p.i])
  prots.grs <- prots.grs[order(names(prots.grs))]
  
  # 5x intersect
  tot.intrsct <- Reduce(GenomicRanges::intersect, prots.grs)
  retl[["int.dom.5x"]] <- tot.intrsct
  
  
  # Subtract 5x intersect from single doms
  prots.grs.wo.intrsct <- lapply(prots.grs, function(lst) GenomicRanges::setdiff(lst, tot.intrsct))
  names(prots.grs.wo.intrsct) <- paste0(names(prots.grs), ".", prot, ".minus.", prot, ".tot.intrsect")
  retl[["leftouts"]] <- prots.grs.wo.intrsct
  
  return(retl)
})

names(totals.and.leftouts) <- prots

# Brains

br.rsem.hp1.cons <- br.rsem
in_tot.intrct <- subsetByOverlaps(tss.gr, totals.and.leftouts$HP1$int.dom.5x)@ranges@NAMES
br.rsem.hp1.cons$HP1.cons <- ifelse(br.rsem$id %in% in_tot.intrct, 1, 0)
out_of.tot.intrct <- subsetByOverlaps(tss.gr,
                                      totals.and.leftouts$HP1$leftouts$BR.HP1.minus.HP1.tot.intrsect)@ranges@NAMES
br.rsem.hp1.cons$HP1.exc.cons <- ifelse(br.rsem$id %in% out_of.tot.intrct, 1, 0)

br.hp1.cons.pval <- format(wilcox.test(br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.cons == 1],
                                br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.exc.cons == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)

pdf("BR.cons.HP1.and.excl.cons.pdf", width = 5)
  boxplot(br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.cons == 1],
          br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.exc.cons == 1],
          outline = F, names = c("In conservative\nHP1 domains", "In HP1 domains\nexcluding conservative"),
          ylab = "TPM", main = "Brains")
  rect(2.25 - 0.3, 80 - 2,
       2.25 + 0.3, 80 + 2)
  
  text(2.25, 80,
       paste0("p = ", br.hp1.cons.pval))
  text(1.3, 36,
       length(br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.cons == 1]))
  text(2.3, 6,
       length(br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.exc.cons == 1]))
dev.off()

boxplot(br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.cons == 1 & br.rsem.hp1.cons$TPM > 0],
        br.rsem.hp1.cons$TPM[br.rsem.hp1.cons$HP1.exc.cons == 1 & br.rsem.hp1.cons$TPM > 0],
        outline = F)

# Fat Bodies

fb.rsem.hp1.cons <- fb.rsem
in_tot.intrct <- subsetByOverlaps(tss.gr, totals.and.leftouts$HP1$int.dom.5x)@ranges@NAMES
fb.rsem.hp1.cons$HP1.cons <- ifelse(fb.rsem$id %in% in_tot.intrct, 1, 0)
out_of.tot.intrct <- subsetByOverlaps(tss.gr,
                                      totals.and.leftouts$HP1$leftouts$FB.HP1.minus.HP1.tot.intrsect)@ranges@NAMES
fb.rsem.hp1.cons$HP1.exc.cons <- ifelse(fb.rsem$id %in% out_of.tot.intrct, 1, 0)

fb.hp1.cons.pval <- format(wilcox.test(fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.cons == 1],
                                fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.exc.cons == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)

pdf("FB.cons.HP1.and.excl.cons.pdf", width = 5)
  boxplot(fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.cons == 1],
          fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.exc.cons == 1],
          outline = F, names = c("In conservative\nHP1 domains", "In HP1 domains\nexcluding conservative"),
          ylab = "TPM", main = "Fat Bodies")
  rect(2.25 - 0.3, 27.5 - 1.1,
       2.25 + 0.3, 27.5 + 1.1)
  
  text(2.25, 27.5,
       paste0("p = ", fb.hp1.cons.pval))
  text(1.3, 12,
       length(fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.cons == 1]))
  text(2.3, 1.2,
       length(fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.exc.cons == 1]))
dev.off()

boxplot(fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.cons == 1 & fb.rsem.hp1.cons$TPM > 0],
        fb.rsem.hp1.cons$TPM[fb.rsem.hp1.cons$HP1.exc.cons == 1 & fb.rsem.hp1.cons$TPM > 0],
        outline = F)

# Kc167

kc.rsem.hp1.cons <- kc.rsem
in_tot.intrct <- subsetByOverlaps(tss.gr, totals.and.leftouts$HP1$int.dom.5x)@ranges@NAMES
kc.rsem.hp1.cons$HP1.cons <- ifelse(kc.rsem$id %in% in_tot.intrct, 1, 0)
out_of.tot.intrct <- subsetByOverlaps(tss.gr,
                                      totals.and.leftouts$HP1$leftouts$Kc167.HP1.minus.HP1.tot.intrsect)@ranges@NAMES
kc.rsem.hp1.cons$HP1.exc.cons <- ifelse(kc.rsem$id %in% out_of.tot.intrct, 1, 0)

kc.hp1.cons.pval <- format(wilcox.test(kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.cons == 1],
                                       kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.exc.cons == 1],
                                       alt = "g")$p.value, digits = 3)


pdf("KC.cons.HP1.and.excl.cons.pdf", width = 5)
  boxplot(kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.cons == 1],
          kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.exc.cons == 1],
          outline = F, names = c("In conservative\nHP1 domains", "In HP1 domains\nexcluding conservative"),
          ylab = "TPM", main = "Kc167")
  rect(2.3 - 0.25, 110 - 2.2,
       2.3 + 0.25, 110 + 2.2)
  
  text(2.3, 110,
       paste0("p = ", kc.hp1.cons.pval))
  text(1.2, 50,
       length(kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.cons == 1]))
  text(2.2, 34,
       length(kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.exc.cons == 1]))
dev.off()

boxplot(kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.cons == 1 & kc.rsem.hp1.cons$TPM > 0],
        kc.rsem.hp1.cons$TPM[kc.rsem.hp1.cons$HP1.exc.cons == 1 & kc.rsem.hp1.cons$TPM > 0],
        outline = F)

# TSS in HP1 and not in Lam or Pc or genes in HP1
# also boxplot of expression

kc.tss.rsem.solo.hp1 <- subsetByOverlaps(tss.gr, domains$Kc167.SoloHP1)@ranges@NAMES
kc.genes.rsem.solo.hp1 <- subsetByOverlaps(genes.gr, domains$Kc167.SoloHP1)@ranges@NAMES
kc.genes.rsem.hp1 <- subsetByOverlaps(genes.gr, domains$Kc167.HP1)@ranges@NAMES

kc.rsem.hp1 <- kc.rsem
kc.rsem.hp1$HP1.genes <- ifelse(kc.rsem$id %in% kc.genes.rsem.hp1, 1, 0)

pdf("KC.genes.HP1.doms.pdf", width = 4)
  boxplot(kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 0],
          kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 1],
          outline = F, names = c("not bound", "bound"),
          ylab = "TPM", main = "Kc167 genes in HP1 domains")
  
  text(1.2, 23,
       length(kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 0]))
  text(2.2, 31.5,
       length(kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 1]))
  text(1.5, 60,
       paste0("p = ",
              format(wilcox.test(kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 0],
                          kc.rsem.hp1$TPM[kc.rsem.hp1$HP1.genes == 1],
                          alt = "l")$p.value, digits = 2)))
dev.off()


# III. Expression in domains and non-domains



# Brains

br.damid.exp$BR.LAM <- factor(br.damid.exp$BR.LAM)
br.damid.exp$BR.HP1 <- factor(br.damid.exp$BR.HP1)
br.damid.exp$BR.PC <- factor(br.damid.exp$BR.PC)

pdf("BR.RNA-seq.expr.in.rel.to.domains.pdf", width = 12)
  op <- par(mfrow = c(1, 3))
  nlam <- tapply(br.damid.exp$TPM, br.damid.exp$BR.LAM, length)
  boxplot(TPM ~ BR.LAM, data = br.damid.exp, outline = F, ylim = c(0, 95),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "LAM")
  text(1.3, 38.6, nlam[1])
  text(2.3, 2, nlam[2])
  text(2.3, 92, 
       paste("p =",
             format(wilcox.test(br.damid.exp$TPM[br.damid.exp$BR.LAM == 0],
                         br.damid.exp$TPM[br.damid.exp$BR.LAM == 1],
                         alt = "g")$p.value, digits = 3, scientific = T)))
  
  nhp1 <- tapply(br.damid.exp$TPM, br.damid.exp$BR.HP1, length)
  boxplot(TPM ~ BR.HP1, data = br.damid.exp, outline = F, ylim = c(0, 95),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "HP1")
  text(1.3, 35.2, nhp1[1])
  text(2.3, 5.6, nhp1[2])
  text(2.3, 92, 
       paste("p =",
             format(wilcox.test(br.damid.exp$TPM[br.damid.exp$BR.HP1 == 0],
                                br.damid.exp$TPM[br.damid.exp$BR.HP1 == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  
  npc <- tapply(br.damid.exp$TPM, br.damid.exp$BR.PC, length)
  boxplot(TPM ~ BR.PC, data = br.damid.exp, outline = F, ylim = c(0, 95),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "PC")
  text(1.3, 28.7, npc[1])
  text(2.3, 5.6, npc[2])
  text(2.3, 92, 
       paste("p =",
             format(wilcox.test(br.damid.exp$TPM[br.damid.exp$BR.PC == 0],
                                br.damid.exp$TPM[br.damid.exp$BR.PC == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  par(op)
dev.off()

# Fat Bodies

fb.damid.exp$FB.LAM <- factor(fb.damid.exp$FB.LAM)
fb.damid.exp$FB.HP1 <- factor(fb.damid.exp$FB.HP1)
fb.damid.exp$FB.PC <- factor(fb.damid.exp$FB.PC)

pdf("FB.RNA-seq.expr.in.rel.to.domains.pdf", width = 12)
  op <- par(mfrow = c(1, 3))
  nlam <- tapply(fb.damid.exp$TPM, fb.damid.exp$FB.LAM, length)
  boxplot(TPM ~ FB.LAM, data = fb.damid.exp, outline = F, ylim = c(0, 32),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "LAM")
  text(1.3, 12.7, nlam[1])
  text(2.3, 0.5, nlam[2])
  text(2.3, 31, 
       paste("p =",
             format(wilcox.test(fb.damid.exp$TPM[fb.damid.exp$FB.LAM == 0],
                                fb.damid.exp$TPM[fb.damid.exp$FB.LAM == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  
  nhp1 <- tapply(fb.damid.exp$TPM, fb.damid.exp$FB.HP1, length)
  boxplot(TPM ~ FB.HP1, data = fb.damid.exp, outline = F, ylim = c(0, 32),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "HP1")
  text(1.3, 9.4, nhp1[1])
  text(2.3, 1.4, nhp1[2])
  text(2.3, 31, 
       paste("p =",
             format(wilcox.test(fb.damid.exp$TPM[fb.damid.exp$FB.HP1 == 0],
                                fb.damid.exp$TPM[fb.damid.exp$FB.HP1 == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  
  npc <- tapply(fb.damid.exp$TPM, fb.damid.exp$FB.PC, length)
  boxplot(TPM ~ FB.PC, data = fb.damid.exp, outline = F, ylim = c(0, 32),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "PC")
  text(1.3, 9.8, npc[1])
  text(2.3, 0.7, npc[2])
  text(2.3, 31, 
       paste("p =",
             format(wilcox.test(fb.damid.exp$TPM[fb.damid.exp$FB.PC == 0],
                                fb.damid.exp$TPM[fb.damid.exp$FB.PC == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  par(op)
dev.off()

# Kc167

kc.damid.exp <- cbind(kc.rsem, sapply(ls(pattern = "Kc167", envir = domains), function(grname){
  gr <- get(grname, envir = domains)
  in_domain <- subsetByOverlaps(tss.gr, gr)@ranges@NAMES
  factor(ifelse(br.rsem$id %in% in_domain, 1, 0))
}))

# write.table(kc.damid.exp, "kc.rsem.domains.csv", quote = F, dec = ",",
#             sep = ";", row.names = F)


pdf("KC.RNA-seq.expr.in.rel.to.domains.pdf", width = 12)
  op <- par(mfrow = c(1, 3))
  nlam <- tapply(kc.damid.exp$TPM, kc.damid.exp$Kc167.LAM, length)
  boxplot(TPM ~ Kc167.LAM, data = kc.damid.exp, outline = F, ylim = c(0, 90),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "LAM")
  text(1.3, 34.7, nlam[1])
  text(2.3, 1.4, nlam[2])
  text(2.3, 86, 
       paste("p =",
             format(wilcox.test(kc.damid.exp$TPM[kc.damid.exp$Kc167.LAM == 0],
                                kc.damid.exp$TPM[kc.damid.exp$Kc167.LAM == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  
  nhp1 <- tapply(kc.damid.exp$TPM, kc.damid.exp$Kc167.HP1, length)
  boxplot(TPM ~ Kc167.HP1, data = kc.damid.exp, outline = F, ylim = c(0, 90),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "HP1")
  text(1.3, 23.7, nhp1[1])
  text(2.3, 31.8, nhp1[2])
  text(2.3, 86, 
       paste("p =",
             format(wilcox.test(kc.damid.exp$TPM[kc.damid.exp$Kc167.HP1 == 0],
                                kc.damid.exp$TPM[kc.damid.exp$Kc167.HP1 == 1],
                                alt = "l")$p.value, digits = 3, scientific = T)))
  
  npc <- tapply(kc.damid.exp$TPM, kc.damid.exp$Kc167.PC, length)
  boxplot(TPM ~ Kc167.PC, data = kc.damid.exp, outline = F, ylim = c(0, 90),
          names = c("not bound", "bound"), ylab = "TPM", xlab = "TSS", main = "PC")
  text(1.3, 29.3, npc[1])
  text(2.3, 1.9, npc[2])
  text(2.3, 86, 
       paste("p =",
             format(wilcox.test(kc.damid.exp$TPM[kc.damid.exp$Kc167.PC == 0],
                                kc.damid.exp$TPM[kc.damid.exp$Kc167.PC == 1],
                                alt = "g")$p.value, digits = 3, scientific = T)))
  par(op)
dev.off()


