#obtain SNPs that exceeds null model of no spatial auto-correlation
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
#setwd("/home/mac/swth/trout")
#setwd("/home/mac/freshwater/combinedBigCreek/")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel")
#setwd("/home/mac/freshwater/combinedBigCreekMal")
setwd(args[1])
load("./out/data.rda")
#load("./out/PermAndObs_10000.rda")
#pob <- PermAndObs_10000
load("./out/PermAndObs_25000.rda")
pob <- PermAndObs_25000
load("./out/SNP_and_Contig_Stats.rda")

zf.cl <- ZF.Chrom.Lengths

autosome.y.max <- 1.0
max.y<-1
max.y.frac = 1
chrom.ht= 1.25

chrom.height <- 1.25  # how much vertical space is given to each chromosome
chrom.space.bottom <- .25  # how much below the fst=0 point each chromosome's frame extends
start.line <- 9
x.start <- 0
max.len <- max(zf.cl$NumBases )
#max(x$Obs$Fst$yhat, x$Perm$Fst[,"100%"]))
# max(pob[[C]]$Obs$Fst$yhat,na.rm=TRUE)
for( C in names(pob)) {
  #print(C)
  y <- SNP.stats[[C]]
  MAF <- (y$Pop1.f1 + y$Pop2.f1) / (y$Pop1.f1 + y$Pop1.f2 + y$Pop2.f1 + y$Pop2.f2)
  #y <-y[MAF>0.05,] #doesn't shorten properly
  fsts <- y$Fst;
  num.snps <- length(fsts)
  pq <- pob[[C]]$Perm$Fst # permutation quantile matrix
  #values for permuted
  kx=pob[[C]]$Obs$Fst$z
  ky.list=list(pob[[C]]$Obs$Fst$yhat, pq[, "99.9%"], pq[, "99.99%"], pq[, "100%"], pq[, "0%"]) # list of kernel density y vectors. [[1]] is observed, [[2]] is 1/1000, [[3]] is 1/10000, [[4]] is max of the upper significance stuff 1 in 25K.
  ky.list[is.na(ky.list)] <- 0
  #ky.list has NA, low SNP density?
  start.line <- 9  #new.vals[[1]] #for me this is non-numeic
 # x.start <- (zf.cl$NumBases[zf.cl$Chromo==C ])+ .07 * max.len  # the .1 here dictates how much space we have after the end of the chromo before the start of the next one
  x.start<-0 # as we aren't plotting, we only want to find the relationship for each chromosome where something is > orange line
  #print(new.vals)
  zero.point <- start.line * chrom.ht
  df<-NULL
  df$Chromosome<-y$Chromosome
  df$GenCoord<-y$GenCoord
  df$POS<-y$POS
  df$idx<-y$CHROM
  df$kx<-kx
  #df$fsts<-fsts # our calculated SNP values...

  #our smoothed fsts are actually
  #lines(x.start + kx, max.y.frac*(ky.list[[1]]/max.y) + zero.point, col="blue", lwd=.7)
  #max.y.frac*(ky.list[[1]]/max.y)
 # df$fsts<-(max.y.frac*(ky.list[[1]]/max.y))
  #df$perm<-(max.y.frac*(ky.list[[4]]/max.y))
  df$perm<-(max.y.frac*(ky.list[[4]]/max.y))
  df$fsts<-(max.y.frac*(ky.list[[1]]/max.y))
  dff<-tbl_df(df)
  SNPS<-dff %>% mutate(diff=fsts-perm) %>% filter(diff > 0)
  #gives us where the blue line > orange line.
  write(C,file="./out/exceededSnps.txt",append=TRUE)
  write.table(SNPS,file="./out/exceededSnps.txt",append=TRUE,sep="\t")
}

#plot(kx, fsts)
#lines(x.start + kx, max.y.frac*(ky.list[[4]]/max.y) + zero.point, col="orange",lwd=.5) # this is at the 1 in 25K level in PermAndObs_25000
#removed zero.point as is omy29
#lines(x.start + kx, max.y.frac*(ky.list[[4]]/max.y), col="orange",lwd=.5) # this is at the 1 in 25K level in PermAndObs_25000
