library(stringr)
library(data.table)
library(dplyr)
library(GenomicRanges)

# Function to shuffle domains
bedTools.shuffle.tss <- function(query, subject, opt.string="-chrom"){
  
  bed.file.1 <- tempfile()
  bed.file.2 <- tempfile()
  
  shuf.1 <- tempfile()
  shuf.2 <- tempfile()
  
  
  options(scipen = 99)
  
  write.table(query, file = bed.file.1, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(subject, file = bed.file.2, quote = F, sep = "\t", col.names = F, row.names = F)
    command = paste("bedtools shuffle -i", bed.file.1,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.1, ";",
                    "bedtools shuffle -i", bed.file.2,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.2)
    # cat(command, "\n")
    try(system(command))
  qu <- import.bed(shuf.1)
  subj <- import.bed(shuf.2)
  q.x.s <- subsetByOverlaps(qu, subj)
  unlink(bed.file.1); unlink(bed.file.2); unlink(shuf.2); unlink(shuf.1)
  return(length(q.x.s))
}



setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/neuro_and_glia_RNAseq/")

# Chromosomes that will be used in analysis
chroms <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

# Read Glia and Neurons rna-seq data
rnaseq <- fread("GSE71104.csv", dec = ",", sep = "\t", stringsAsFactors = F)

# read list of ubiquitous genes
ub.genes <- read.csv("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/RNAseq_vs_DamID/Ubiq_genes/ubiqex.genes.csv", stringsAsFactors = F)

# fix obnoxious formatting of chromosome coordinates 
coords <- as_data_frame(str_split_fixed(rnaseq$`chrloc (dm3)`, ":|-", 3), stringAsFactors = F) %>%
  setNames(c("chr", "start", "end")) %>% 
  mutate_each_(funs(as.integer), 2:3) %>% 
  mutate(TSS = ifelse(rnaseq$strand == 1, start, end))

# Subsetting of genes, expressing in neurons, chromosome filtering
neurnaseq <- data.frame(
  FBGN = rnaseq$ensembl,
  name = rnaseq$symbol,
  coords,
  strand = rnaseq$strand,
  TPM = (rnaseq$Neuron_rep1 + rnaseq$Neuron_rep2) / 2,
  biotype <- rnaseq$biotype
) %>% mutate(TSS = ifelse(strand == 1, start, end)) %>% 
               filter(biotype == "protein_coding", TPM > 1,
                      !(FBGN %in% ub.genes$FlyBase.ID)) %>%
  filter(chr %in% chroms)

# Make GRanges object with TSS of above genes
neur.tss.gr <- GRanges(
  seqnames = Rle(neurnaseq$chr),
  ranges = IRanges(
    start = neurnaseq$TSS,
    width = 1, 
    names = neurnaseq$FBGN
  )
)

# Subsetting of genes, expressing in glia, chromosome filtering
glirnaseq <- data.frame(
  FBGN = rnaseq$ensembl,
  name = rnaseq$symbol,
  coords,
  strand = rnaseq$strand,
  TPM = (rnaseq$Glia_rep1 + rnaseq$Glia_rep2) / 2,
  biotype <- rnaseq$biotype
) %>% mutate(TSS = ifelse(strand == 1, as.integer(start), as.integer(end))) %>% 
  filter(biotype == "protein_coding", TPM > 1, !(FBGN %in% ub.genes$FlyBase.ID),
         chr %in% chroms)


# Make GRanges object with TSS of above genes
glia.tss.gr <- GRanges(
  seqnames = Rle(glirnaseq$chr),
  ranges = IRanges(
    start = glirnaseq$TSS,
    width = 1, 
    names = glirnaseq$FBGN
  )
)

# folder, contatining .beds of domains, identified by HMM
dom.dir <- paste0("~/IMG/Projects/",
                  "HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/",
                  "DamID-seq.HP1.PC.Lam.WBr.Nrn.Glia.Fb/final_variant/",
                  "BioHMM2.qn.full.PC.HMM3/")
# Create new environment and load domains there
domains <- new.env()
for (i in dir(dom.dir, pattern = "(NRN|Glia).*bed")){
  bed <- fread(file.path(dom.dir, i), skip = 1) %>% 
    setNames(c("chr", "start", "end"))
  assign(
    sub("\\.domains\\.bed", "", i),
    GRanges(
      seqnames = Rle(bed$chr),
      ranges = IRanges(
        start = bed$start,
        end = bed$end
      )
    ),
    envir = domains
  )
}
rm(i, bed)
# Subset neurone-specific genes by their presence/absence in domains
# and add corresponding column to rnaseq dataset
# Then plot boxplots and 
nrn.x.lam <- subsetByOverlaps(neur.tss.gr, domains$NRN.LAM)
neurnaseq$LAM <- ifelse(neurnaseq$FBGN %in% names(nrn.x.lam), 1, 0)

boxplot(TPM ~ LAM, neurnaseq, outline = F)
wilcox.test(TPM ~ LAM, neurnaseq, alt = "g")


nrn.x.hp1 <- subsetByOverlaps(neur.tss.gr, domains$NRN.HP1)
neurnaseq$HP1 <- ifelse(neurnaseq$FBGN %in% names(nrn.x.hp1), 1, 0)

boxplot(TPM ~ HP1, neurnaseq, outline = F)
wilcox.test(TPM ~ HP1, neurnaseq, alt = "g")

nrn.x.pc <- subsetByOverlaps(neur.tss.gr, domains$NRN.PC)
neurnaseq$PC <- ifelse(neurnaseq$FBGN %in% names(nrn.x.pc), 1, 0)

boxplot(TPM ~ PC, neurnaseq, outline = F)
wilcox.test(TPM ~ PC, neurnaseq, alt = "g")

sapply(neurnaseq[, 10:12], function(x) sum(x)/length(x))
sapply(neurnaseq[, 10:12], function(x) (length(x) - sum(x))/length(x))


glia.x.lam <- subsetByOverlaps(glia.tss.gr, domains$Glia.LAM)
glirnaseq$LAM <- ifelse(glirnaseq$FBGN %in% names(glia.x.lam), 1, 0)

boxplot(TPM ~ LAM, glirnaseq, outline = F)
wilcox.test(TPM ~ LAM, glirnaseq, alt = "g")


glia.x.hp1 <- subsetByOverlaps(glia.tss.gr, domains$Glia.HP1)
glirnaseq$HP1 <- ifelse(glirnaseq$FBGN %in% names(glia.x.hp1), 1, 0)

boxplot(TPM ~ HP1, glirnaseq, outline = F)
wilcox.test(TPM ~ HP1, glirnaseq, alt = "g")

glia.x.pc <- subsetByOverlaps(glia.tss.gr, domains$Glia.PC)
glirnaseq$PC <- ifelse(glirnaseq$FBGN %in% names(nrn.x.pc), 1, 0)

boxplot(TPM ~ PC, glirnaseq, outline = F)
wilcox.test(TPM ~ PC, glirnaseq, alt = "g")

sapply(glirnaseq[, 10:12], function(x) sum(x)/length(x))
sapply(glirnaseq[, 10:12], function(x) (length(x) - sum(x))/length(x))

glirnaseq$TPM.cut <- cut(glirnaseq$TPM, quantile(glirnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
                         labels = c("weak", "mid", "high"), include.lowest = T)


nindoms <- sapply(glirnaseq  %>% select(10:12),  function(x) sum(x == 0))



interdom.ratios <- as.data.frame(sapply(levels(glirnaseq$TPM.cut), function(i){
  sapply(glirnaseq %>% filter(TPM.cut == i) %>% select(10:12),  function(x) sum(x == 0)/length(x))
}))

interdom.ratios$n <- nindoms

write.table(interdom.ratios, "specific.genes.in.interdoms.by.exp.glia.csv", sep = "\t", dec = ',')

neurnaseq$TPM.cut <- cut(neurnaseq$TPM, quantile(neurnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
                         labels = c("weak", "mid", "high"), include.lowest = T)


nindoms <- sapply(neurnaseq[, 10:12],  function(x) sum(x == 0))



interdom.ratios <- as.data.frame(sapply(levels(neurnaseq$TPM.cut), function(i){
  sapply(neurnaseq[neurnaseq$TPM.cut == i, 10:12],  function(x) sum(x == 0)/length(x))
}))

interdom.ratios$n <- nindoms

write.table(interdom.ratios, "specific.genes.in.interdoms.by.exp.neur.csv", sep = "\t", dec = ',', quote = F)
sapply(neurnaseq[neurnaseq$TPM.cut == "weak", 10:12],  function(x) length(x))

write.table(neurnaseq, "neurnaseq_specific_genes_lam_hp1_pc_enrich.csv", sep = "\t", quote = F, row.names = F)
write.table(glirnaseq, "glirnaseq_specific_genes_lam_hp1_pc_enrich.csv", sep = "\t", quote = F, row.names = F)


overlap.pval <- function(rnaseq, doms){
  tss.bed <- glirnaseq[, c(3,6)] %>% mutate(end = TSS + 1) %>% arrange(chr, TSS)
  tss.gr <- GRanges(seqnames = Rle(tss.bed$chr), ranges = IRanges(start = tss.bed$TSS, width = 1))
  gaps.doms <- gaps(doms)
  gaps.bed <- data.frame(chr = seqnames(gaps.doms),
                         start = start(gaps.doms),
                         end = end(gaps.doms))
  et <- length(subsetByOverlaps(tss.gr, gaps.doms))
  shuff <- sapply(1:10000, function(i) bedTools.shuffle.tss(tss.bed, gaps.bed))
  return(sum(shuff >= et) / 10000)
}

glia.pval <- sapply(names(domains)[grepl("Glia", names(domains))],
       function(nom){
         overlap.pval(glirnaseq, domains[[nom]])
       })
overlap.pval(glirnaseq, domains$Glia.PC)


tss.bed <- glirnaseq[, c(3,6)] %>% mutate(end = TSS + 1) %>% arrange(chr, TSS)
tss.gr <- GRanges(seqnames = Rle(tss.bed$chr), ranges = IRanges(start = tss.bed$TSS, width = 1))
gaps.hp1 <- gaps(domains$Glia.HP1)
length(subsetByOverlaps(tss.gr, gaps.hp1))

gaps.hp1.bed <- data.frame(chr = seqnames(gaps.hp1),
                           start = start(gaps.hp1),
                           end = end(gaps.hp1))
bedTools.shuffle.tss(tss.bed, gaps.hp1.bed)

nrn.pval <- sapply(names(domains)[grepl("NRN", names(domains))],
                    function(nom){
                      overlap.pval(neurnaseq, domains[[nom]])
                    })


plot(density(as.numeric(neurnaseq$end) - as.numeric(neurnaseq$start), bw = 100), xlim = c(0, 15000))

plot(density(as.numeric(glirnaseq$end) - as.numeric(glirnaseq$start)), xlim = c(0, 15000))


# Count tissue-specific genes, which associated with negative damid values



