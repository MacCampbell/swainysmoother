#! /usr/bin/Rscript

#the goal is to reduce windows we aren't interested in, >150 bp and < 1000 bp

args = commandArgs(trailingOnly=TRUE)

file<-args[1]
outfile<-args[2]

data<-read.table(file)
data$ranges<-data$V3-data$V2

data2<-subset(data, ranges>150 & ranges < 1000)

write.table(data2, file=paste(outfile), sep=" ")


