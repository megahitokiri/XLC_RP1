#source("https://bioconductor.org/biocLite.R")
#biocLite("BioSeqClass")
  #biocLite("sangerseqR")

library(knitr, quietly=TRUE)
library(sangerseqR, quietly=TRUE)
library(Biostrings, quietly=TRUE)
opts_chunk$set(tidy=TRUE)

Sanger_data <- "C:/Users/megah/OneDrive/Espana Paper/Gencove_Retinitis/Sanger/RP1_JLG.ab1"

RP1_JLG <- read.abif(Sanger_data)

str(RP1_JLG, list.len=20)

Sanger_JLG <- sangerseq(RP1_JLG)
#Seq1 <- primarySeq(Sanger_JLG)

primarySeq(Sanger_JLG, string = TRUE)

chromatogram(Sanger_JLG, width = 200, height = 2, trim5 = 50, trim3 = 100, showcalls = "both",
             filename = "chromatogram.pdf")

chromatogram(Sanger_JLG, width = 100, height = 2, trim5 = 0, trim3 = 100, showcalls = "primary",
             filename = "chromatogram2.pdf")
  