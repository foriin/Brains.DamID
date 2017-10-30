library(data.table)
library(dplyr)
library(tidyr)

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/HP1.kd.MA.Kc/")
chip.key <- fread("affymetrix_key_GPL1322-26772.txt")%>%
  select(1,2,10,11) %>% 
  rename(FBgn = CLONE_ID_LIST)

gene.coor <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf")

gene.coor <- gene.coor %>% 
  select(1,2,4,5,7,9) %>% 
  setNames(c("chr", "source", "start", "end", "strand", "mess")) %>% 
  separate(mess, c("FBgn", "trash", "gene_name"), ";") %>% 
  mutate(FBgn = sub('.*"(.*)"', "\\1", FBgn), gene_name = sub('.*"(.*)"', "\\1", gene_name)) %>% 
  select(-trash)

chip.key.r5 <- merge(chip.key, gene.coor, by = "FBgn")

gse18092.r5 <- merge(chip.key.r5, ma_data2, by = "ID") %>% select(1:5, 7:9, 17:19) %>% 
  mutate(TSS = ifelse(strand == "+", start, end))

write.table(gse18092.r5, "GSE18092.r5.gene.coord.tss.csv", row.names = F, sep = '\t')
