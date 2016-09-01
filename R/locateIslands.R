# this is a script that goes through the PermAndObs25000 object and
# delineates the start and stop positions of islands on each chromosome
# and stores those in a list of data frames
#Mac in my case it is  PermAndObs_10000.rda

# Set to your own working directory:
#setwd("/Users/eriq/Documents/work/assist/kristenruegg/SWTH_genome_analysis_for_paper")
setwd("/home/mac/freshwater/combinedBigCreek")
# load the data we need:
#load("inputs_and_scripts/swth_snps_mapped_to_ZF.rda")
load("data.rda")
rownames(ZF.Chrom.Lengths) <- ZF.Chrom.Lengths$Chromo  # give these some rownames

#load("inputs_and_scripts/PermAndObs_25000_all.rda")
load("./swainy/PermAndObs_10000.rda")
#for convenience, renaming any object to this

#load("inputs_and_scripts/SNP_and_Contig_Stats.rda")  # We will want this to say which SNPs are in islands, etc
load("SNP_and_Contig_Stats.rda") #output from BasicPopGenStatistics.R

#PermAndObs_25000[["Z"]]$Obs$Fst$yhat[is.na(PermAndObs_25000[["Z"]]$Obs$Fst$yhat)] <- 1.0
#PermAndObs_25000[["Z"]]$Perm$Fst[is.na(PermAndObs_25000[["Z"]]$Perm$Fst)] <- 1.0
#Mac I don't have a sex chromosome

## here we have to do a little data cleaning: there was one stretch on the z chromosome where the
## smoother returned NaN, most likely due to data sparsity.  It is right next to an observed smooth
## block of 1.0.  In the Perm results there are some NAs too, but the "1 in 25K" value is 1.0 anyway
## nearby.  So, I will just replace those NAs with 1.0 in the Fst results (it is about 6 smoothed points
## out of over 4100


# now define a function.  We give it:
#  x, 	the x-values at which smoothed predicted values are available,
#  y, 	the predicted values at each x
#  yc, 	the "critical values" of y at each x.  That is, if y>yc then that
#				spot is in an island
#
#
#  AND, this function returns a vector of start and stop points in the
# x-dimension.  The distance between successive x values should be constant, call it d
# Islands are deemed to start at the greater of the first (x value where y>yc)-0.5d or 0
# and they start 0.5*d above the last one (or said a different way, below the first FALSE
# after a TRUE, making sure it does not exceed the chromo length
#
StartAndEndPoints <- function(x, y, yc, max.x) {
  #Mac added this to see if there are any values

  d <- x[2]-x[1]

  st <- y>yc  # vector of states (TRUE if y>yc, FALSE otherwise)

  # now cast each spot as a pair, i.e. going from a state to another state.
  # we define the 0 position from as FALSE and the last position to as FALSE
  from <- c(F,st)
  to <- c(st,F)

  starts <- c(0,x)[from!=to & to==TRUE]
  ends <- c(x,max.x)[from!=to & to==FALSE]

  # now make the starts 1/2 * d less
  starts <- starts - d/2
  starts[starts<0] <- 0

  # and make the ends d/2 greater
  ends <- ends + d/2
  ends[ends>max.x] <- max.x

  data.frame(Starts=starts, Ends=ends)
}


### now we just need to apply these to all the chromosomes:
thresh <- "100%"  # here we are pulling it out at the 1 in 25,000 level
#going for the threshold of 100% and changing names
#IslandStartAndEnds <- lapply(names(PermAndObs_10000), function(h) {
#  print(H)
#  StartAndEndPoints(
#    x=PermAndObs_10000[[h]]$Obs$Fst$z,
#    y=PermAndObs_10000[[h]]$Obs$Fst$yhat,
#    yc=PermAndObs_10000[[h]]$Perm$Fst[,thresh],
#    max.x=ZF.Chrom.Lengths[h,"NumBases"]
#  )
#}
#)

#From Mac
IslandStartAndEnds<-lapply(names(PermAndObs_10000), function(h) {
  StartAndEndPoints(
  x=PermAndObs_10000[[h]]$obs$Fst$z,
  y=PermAndObs_10000[[h]]$Obs$Fst$yhat,
  yc=PermAndObs_10000[[h]]$Perm$Fst[,thresh],
  max.x=ZF.Chrom.Lengths[h,"NumBases"]
  )
}
)

names(IslandStartAndEnds) <- names(PermAndObs_10000)


### and it turns out I wouldn't mind having the same things for
### the 99.9% and 99.99% cutoffs too (1 in 1,000 and 1 in 10,000)
IslandStartEndsList <- list()
#IslandStartEndsList$"OneIn25K" <- IslandStartAndEnds
#Mac, 1 in 10K for this set
IslandStartEndsList$"OneIn10K" <- IslandStartAndEnds


# do it for 1 in 10,000
#thresh <- "99.99%"
#IslandStartEndsList$"OneIn10K" <- lapply(names(PermAndObs_25000), function(h) {
 # print(h)
#  StartAndEndPoints(
    #x=PermAndObs_25000[[h]]$Obs$Fst$z,
    #y=PermAndObs_25000[[h]]$Obs$Fst$yhat,
   # yc=PermAndObs_25000[[h]]$Perm$Fst[,thresh],
  #  max.x=ZF.Chrom.Lengths[h,"NumBases"]
 # )
#}
#)
#names(IslandStartEndsList$"OneIn10K") <- names(PermAndObs_25000)


# and finally for 1 in 1,000
#thresh <- "99.9%"
#IslandStartEndsList$"OneIn1K" <- lapply(names(PermAndObs_25000), function(h) {
  #print(h)
  #StartAndEndPoints(
  #  x=PermAndObs_25000[[h]]$Obs$Fst$z,
  #  y=PermAndObs_25000[[h]]$Obs$Fst$yhat,
 #   yc=PermAndObs_25000[[h]]$Perm$Fst[,thresh],
#    max.x=ZF.Chrom.Lengths[h,"NumBases"]
# )
#}
#)
#names(IslandStartEndsList$"OneIn1K") <- names(PermAndObs_25000)






### and now we can get the proportion of each chromosome that is in islands
PropInIslands <- lapply(
  names(IslandStartAndEnds),
  function(x) {
    sum(IslandStartAndEnds[[x]]$Ends-IslandStartAndEnds[[x]]$Starts) / ZF.Chrom.Lengths[x,"NumBases"]
  }
)
names(PropInIslands) <- names(IslandStartAndEnds)




### and now we can say whether every SNP is in an island or not
SNP.stats.w.islands <- lapply(names(IslandStartAndEnds), function(h) {
  InIsland <- sapply(SNP.stats[[h]]$GenCoord, function(x) as.logical(sum(IslandStartAndEnds[[h]]$Starts<x & IslandStartAndEnds[[h]]$Ends>x)) )
  data.frame(SNP.stats[[h]], InIsland)
}
)
names(SNP.stats.w.islands) <- names(SNP.stats)

### and we can do the same thing for contigs too.  I am going to say that a contig is in an island if its
### average location is in the island
contig.stats.w.islands <- lapply(names(IslandStartAndEnds), function(h) {
  InIsland <- sapply(contig.stats[[h]]$GenCoord, function(x) as.logical(sum(IslandStartAndEnds[[h]]$Starts<x & IslandStartAndEnds[[h]]$Ends>x)) )
  data.frame(contig.stats[[h]], InIsland)
}
)
names(contig.stats.w.islands) <- names(IslandStartAndEnds)



### and now we can compute mean and sd of statistics for things in islands and not in islands.
### I will want to do this for the autosomes and Z chrom separately. Will have to fix this when I have the Z
se.mean <- function(x) { sqrt(var(x)/length(x)) }
noz <- SNP.stats.w.islands[names(SNP.stats.w.islands)!="Z"] # list of snp stats with no Z chrom
wz <- SNP.stats.w.islands[names(SNP.stats.w.islands)=="Z"] # list of snp stats that is just the Z chrom
noz <- do.call(rbind, noz)  # make it a data frame
wz <- do.call(rbind, wz)



# now do the calcs for Fst, the only SNP-specific measure:
means <- rbind(c(t(aggregate(Fst ~ InIsland, noz, mean))[2,],  t(aggregate(Fst ~ InIsland, wz, mean))[2,]))
ses <- rbind(c(t(aggregate(Fst ~ InIsland, noz, se.mean))[2,],  t(aggregate(Fst ~ InIsland, wz, se.mean))[2,]))
rownames(means)<-c("Fst")
rownames(ses)<-c("Fst")
colnames(means) <- c("autosom.bg", "autosom.in.islands", "z.chrom.bg", "z.chrom.in.islands")
colnames(ses) <- colnames(means)


### and let us do the same for contigs
ccnoz <- contig.stats.w.islands[names(contig.stats.w.islands)!="Z"]
ccwz <- contig.stats.w.islands[names(contig.stats.w.islands)=="Z"]
ccnoz <- do.call(rbind,ccnoz)
ccwz <- do.call(rbind, ccwz)

means <- rbind(means,
               cbind(t(aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ InIsland, ccnoz, mean)), t(aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ InIsland, ccwz, mean)))
)

ses <- rbind(ses,
             cbind(t(aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ InIsland, ccnoz, se.mean)), t(aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ InIsland, ccwz, se.mean)))
)


# finally, I think it would be good to get the total proportion of genome in and
# out of blocks in autosomes and also on the Z chromo
len <- ZF.Chrom.Lengths$NumBases
names(len) <- ZF.Chrom.Lengths$Chromo
autos<-names(PropInIslands)[names(PropInIslands)!="Z"]

pp.auto <- weighted.mean(unlist(PropInIslands[autos]), len[autos])  # proportion of autosomes in islands
pp.z <- PropInIslands[["Z"]]
pp.g <- rbind(c(pp.auto, 1-pp.auto, pp.z, 1-pp.z))
rownames(pp.g) <- "Ppn.Of.Genome"
colnames(pp.g) <- colnames(means)


# here are the rows we want to print out of those:
rrs <- c("df", "Fst", "Pi.2", "Pi.1", "dxy")

# here we can print it out:
format(means[rrs,], digits=2)
format(ses[rrs,], digits=2)


# catenate the list into a big data frame:
SNP.stats.with.islands <- do.call(rbind, args=SNP.stats.w.islands)



### and here we can spit out the data frame, both as an Rda and as a text file
## so that we can run snps_falling_in_genes.sh on it.
save(SNP.stats.with.islands, file="SNP.stats.with.islands.rda")
write.table(SNP.stats.with.islands, col.names=T, row.names=F, sep=" ", quote=F, file="SNP.stats.with.islands.txt")

# and save the IslandStartAndEnds as a list to put in the multi-chrom figure later
save(IslandStartAndEnds, file="IslandStartAndEnds.rda")
save(IslandStartEndsList, file="IslandStartEndsList.rda")  # save the list too.  We will use this for the supplemental figures.


# and as a final hurrah let us get the min, max, and average lengths of these islands
ai <- do.call(rbind, IslandStartAndEnds[names(IslandStartAndEnds)!="Z"])  # ai ==> autosome islands
ai$Length <- ai$Ends - ai$Starts
c(min(ai$Length), max(ai$Length), mean(ai$Length))

zi <- do.call(rbind, IslandStartAndEnds[names(IslandStartAndEnds)=="Z"])  # zi ==> Z-chromosome islands
zi$Length <- zi$Ends - zi$Starts
c(min(zi$Length), max(zi$Length), mean(zi$Length))




# just playing around here to see if I can make a plot of the sort that
# kristen wanted, but this is way lame....
#par(mar=c(.5, 2, .5, .5))
#par(mfrow=c(5,1))
#for(r in rrs) {
#plot(1:4, means[r,], xaxt="n", bty="l", xlab="")
#segments(1:4, means[r,]-2.0*ses[r,], 1:4, means[r,]+2.0*ses[r,], )
#}
