#reads in "convertedData.txt" from swainyPrepper.R and makes it into an .rda with naming compatible with swainysmoother
#needs a text file with chromosome lengths, here "RT.Chrom.Lengths"
#setwd("/home/mac/swth/trout/dat")
#setwd("/home/mac/freshwater/combinedBigCreek")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel")
setwd("/home/mac/freshwater/combinedBigCreekMal")
library(dplyr)

table<-rename(tbl_df(read.table("./convertedData.txt")), POS=site, Chromosome=chrom, GenCoord=left, CHROM=contig_idx)
ZF.Chrom.Lengths<-read.table("./RT.Chrom.Lengths", header=TRUE)
swth.SNPs<-select(table,Chromosome,GenCoord, CHROM, length, POS, Pop1.f1, Pop1.f2, Pop2.f1, Pop2.f2)

save(swth.SNPs, ZF.Chrom.Lengths, file="./data.rda")
