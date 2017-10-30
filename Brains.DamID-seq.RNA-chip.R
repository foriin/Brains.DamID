

rm(list = ls())

setwd("~/IMG/Projects/HP1.Lamin.Polycomb.DNA.contacts.Effect.on.expression/rnaseq_ns/ExpressionData/")
desalvo <- fread("DeSalvo 2014 Front Neurosci (Suppl)/Transcriptome of glia, surface glia, neurons and brains from imago_deSalvo_2014.csv", dec = ",")
br.rna.chip <- desalvo[, c(7,8,9,3,4,5,11:16,10)]
names(br.rna.chip)[13] <- "strand"
