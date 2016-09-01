# this is a script that uses the PermAndObsXXXXobject to
# plot all Fst smooth results for all chromosomes on a
# single page.

library(plyr)
#pdf("test.pdf")

# Set to your own working directory:
#setwd("/Users/eriq/Documents/work/assist/kristenruegg/SWTH_genome_analysis_for_paper")
#setwd("/home/mac/freshwater/combinedBigCreek")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel")
setwd("/home/mac/freshwater/combinedBigCreekMal")
# load the data we need:
#load("inputs_and_scripts/swth_snps_mapped_to_ZF.rda")
load("data.rda")

rownames(ZF.Chrom.Lengths) <- ZF.Chrom.Lengths$Chromo  # give these some rownames

#load("inputs_and_scripts/PermAndObs_25000_all.rda")
load("./swainy/PermAndObs_10000.rda")
#for variable names
pob <- PermAndObs_10000

# and get the snp Fst stats etc:
#load("inputs_and_scripts/SNP_and_Contig_Stats.rda")
load("SNP_and_Contig_Stats.rda")

# load up the locations where island start and end
#load("IslandStartAndEnds.rda")
#Doesn't apply...




################ FUNCTION DEFINITIONS  ############
plot.chromo.sparkline <- function(
								chrom, # string name of the chromosome (like "1" or "2")
								chrom.len, # number of bases of this chromosome
								pos, # the genome coordinates (position on chrom) of the SNP
								fsts, # the fst values of the SNPs
								chrom.ht,  # vertical space given to each chromosome.  start.line * chrom.ht gives the fst=0 point for the current line
								x.span, # the number of bases represented as the max distance in the horiz direction (basically length of longest chromosome)
								start.line=16, # the line up from the bottom that we will try to print this on
								x.start=0,  # x value at which we want this chrom's  plot to be placed at.  Note, that we must account for extra space needed to the left of it for annotations outside of this function
								kx=NULL, # the kernel density x vector
								ky.list=NULL, # list of kernel density y vectors. [[1]] is observed,
											#[[2]] is upper 1/100, [[3]] is upper 1/1000, [[4]] is upper
											# and [[5]] is lower 1/1000 quantiles of smooths of permuted fst values
								y1.kick=7.5e5, # how far to the left to kick over the left y axis
								plot.genes=F,  # should be we plot some gene positions (given in gene.snps.to.plot)
								max.y = 1,  # set this to the max smoothed or observed value in all the chromosomes to which this is comparable (bascially autosomes or Z chrom)
								max.y.frac = 1,  # how how up on the Fst scale should that max.y value come as a fraction of itself?
								StartEnds = data.frame(NULL)  # should have a Starts column and an Ends column. Shows start and end points of islands in that chromo.
								) {


	# first check to see if we can fit this thing on the current line.  If so, we do it,
	# otherwise we decrement start.line by one (and we return it, too, so that we
	# can follow up with the next plot, etc.
	if( (x.start + chrom.len) > x.span ) {
		start.line <- start.line - 1
		x.start <- 0
	}
	zero.point <- start.line * chrom.ht
	x.end=x.start+chrom.len



	segments(x.start, zero.point + 1, x.start + chrom.len, zero.point + 1, lwd=.2, col="yellow") # draw the upper yellow line before the gray dots
	segments(x.start-y1.kick, zero.point, x.start-y1.kick, zero.point+1, lwd=.3)  # this is like a y-axis. We put it a bit to the left of the start of the data
  #Mac since I had so few I made it black
#	points(x.start + pos, zero.point + fsts, cex=.05, pch=19, col="lightgray")
	points(x.start + pos, zero.point + fsts, cex=.1, pch=19, col="black")
	segments(x.start-y1.kick, zero.point, x.start + chrom.len + y1.kick, zero.point, lwd=1, col="black") # This is like an x-axis
	text(x.start, zero.point + .3 * chrom.ht, labels=C, pos=2, cex=.75)

	# here are tick marks on the y-axis for gray points
	segments(x.start-y1.kick, zero.point + (1:4)/4, x.start-1.5*y1.kick, zero.point + (1:4)/4, lwd=.3)
	text(x=x.start+0.2*y1.kick, y= zero.point + 1, pos=2, cex=0.5, labels="1.0")

	# here is the right hand axis:
	segments(x.start + chrom.len + y1.kick, zero.point, x.start + chrom.len + y1.kick, zero.point+1, lwd=.3)

	# here are the right hand axis tick marks for non-z-chromos:
	if(C=="Z") {
		segments(x.start + chrom.len + y1.kick, zero.point + (1:4)/4, x.start + chrom.len + 1.5*y1.kick, zero.point + (1:4)/4, lwd=.3)
		text(x=x.start + chrom.len - .05*y1.kick, y=zero.point + 1, labels="1.0", pos=4, cex=0.5)  # and the 1.0 text label for the top tick marks
	} else {
		ybits <- zero.point + max.y.frac*(1/4)/max.y
		segments(x.start + chrom.len + y1.kick, ybits , x.start + chrom.len + 1.5*y1.kick, ybits, lwd=.3)
		text(x=x.start + chrom.len - .05*y1.kick, y=ybits, labels="0.25", pos=4, cex=0.5)  # 0.25 text label on the tick mark
	}

	# here is an experimental blue blob that shows where the max blue spot is on the 0 to 1 Fst scale:
	#segments(x.start + chrom.len + 2.2*y1.kick, zero.point , x.start + chrom.len + 2.2*y1.kick, zero.point+max(ky.list[[1]], na.rm=T), lwd=2, lend="butt", col="blue")

	# now start drawing the smoothed lines
	if( !is.null(kx) && !is.null(ky.list)) {
		# plot the significance lines at the 1/1000 level
	  #MAC editing so that NaN -> 1, I don't like it as much and changed it back
		lines(x.start + kx, max.y.frac*(ky.list[[4]]/max.y) + zero.point, col="orange",lwd=.5) # this is at the 1 in 25K level in PermAndObs_25000
	#  permLine<-max.y.frac*(ky.list[[4]]/max.y) + zero.point
	 # permLine[is.na(permLine)]<-1
	  #lines(x.start + kx, permLine, col="orange",lwd=.5)
	#		lines(x.start + kx, ky.list[[5]] + zero.point, col="black",lwd=.5)

		# now draw the observed smooth over that
		lines(x.start + kx, max.y.frac*(ky.list[[1]]/max.y) + zero.point, col="blue", lwd=.7)

		# now point out where the observed smoothed values exceed the significance thresholds
		#ks.above.thresh(kx, ky.list[[1]], ky.list[[2]], x.add=x.start, plot.ht=zero.point + 1.2, col="green", lwd=1.3)
		#above.K <- ks.above.thresh(kx, ky.list[[1]], ky.list[[3]], x.add=x.start, plot.ht=zero.point + 1.35, col="orange", lwd=1.3)
		#ks.above.thresh(kx, ky.list[[1]], ky.list[[4]], x.add=x.start, plot.ht=zero.point + 1.5, col="hotpink1", lwd=1.3)

		# now print, the fraction of genome that is in areas with smoothed fst above the 1/1000 quantile
		#if(!is.null(above.K)) {
	#		fract.above <- sum(above.K[,2] - above.K[,1]) / chrom.len
		#}
		#else {
	#		fract.above <- 0
		#}
		#text(x.end - .014 * x.span, zero.point + .6 * chrom.ht, labels=paste( round(fract.above*100), "%", sep=""), pos=4, cex=.5, col="orange")
	}
	# MAC, we don't have these
	# here we plot the gene positions:
	#if(plot.genes==TRUE) {
	#	gg <- gene.snps.to.plot[ gene.snps.to.plot$Chromosome==C & gene.snps.to.plot$GeneFunctionCategory=="Migration",]
	#	if(nrow(gg)>0) {
	#		points(x.start + gg$GenCoord, zero.point + gg$Fst, #rep(1, length(gg$GenCoord)),
	#			pch=17,
	#			col="deeppink2",
	#			cex=.85
	#		)
	#	}
	#}
	# Mac, we don't have islands
	# here we put the island bars...
	#BH <- 1.17  # for the height of the island bars...
	#if(nrow(StartEnds)>0) {
	#	segments(x.start + StartEnds$Starts, zero.point + BH,  x.start + StartEnds$Ends, zero.point + BH, #lend="butt", lwd=3.5, col="blue")
#	}

	return(list(start.line, x.end))
}
####################### END FUNCTION DEFINITIONS  ###################

#MAC don't have these yet so commenting out
## read in the gene information:
#gene.df <- read.table("SNP_in_Genes_of_Interest.txt", header=T, sep="\t")
# now figure out the max fsts in each of those genes of interest:
#gene.max.fst.df <- ddply(gene.df, "Ensembl.Gene.ID", function(x) x[which.max(x$Fst),] )
#gene.snps.to.plot <- gene.max.fst.df[, c("Chromosome", "GenCoord", "GeneFunctionCategory", "Fst")]
#gene.snps.to.plot$GeneFunctionCategory <- factor(gene.snps.to.plot$GeneFunctionCategory, levels=c("Migration", "Song", "Color"))

# some nice short variable names:
zf.cl <- ZF.Chrom.Lengths



# some stuff for a MAF cutoff for which SNPs to show (do later)

# set some plot area parameters, etc:
chrom.height <- 1.25  # how much vertical space is given to each chromosome
chrom.space.bottom <- .25  # how much below the fst=0 point each chromosome's frame extends
start.line <- 9
x.start <- 0
max.len <- max(zf.cl$NumBases )

par(mar=c(1,1,1,1))
plot(c(0,max.len), c(-.25, (start.line+1) * chrom.height), type="n", axes=F, xlab="", ylab="")
#abline(h=((0:11) * chrom.height) - chrom.space.bottom, lwd=.5)  # black lines to separate rows

# catch the autosomes to compute some values
autosomes <- names(pob)[names(pob)!="Z"]
# the following are for setting the scale of the y-axis for the smooths
#autosome.y.max <- max(sapply(pob[autosomes], function(x) max(x$Obs$Fst$yhat, x$Perm$Fst[,"100%"]))) #MAC NA's in here ruining it all
#autosome.y.max <-0.75 #heh
autosome.y.max <- max(sapply(pob[autosomes], function(x) max(x$Obs$Fst$yhat, x$Perm$Fst[,"100%"],na.rm=TRUE))) #MAC NA's in here ruining it all

#z.chrom.y.max <- sapply(pob["Z"], function(x) max(x$Obs$Fst$yhat, x$Perm$Fst[,"100%"], na.rm=T))  #MAC no z.chrom for us

# now we cycle over chromosomes and plot them out. The plot function
# returns the line that it drew the plot on etc

for( C in names(pob)) { #autosomes) {
#for ( C in c("omy01")) {
#for ( C in c("omy11","omy12","omy13","omy14","omy15","omy16","omy17","omy18","omy19","omy20")) {
#for ( C in c("omy21","omy22","omy23","omy24","omy25","omy26","omy27","omy28","omy29")) {
	y <- SNP.stats[[C]]  # I could add a filter on MAFs here if desired
	MAF <- (y$Pop1.f1 + y$Pop2.f1) / (y$Pop1.f1 + y$Pop1.f2 + y$Pop2.f1 + y$Pop2.f2)
	#MAF[MAF>0.5] <- 1.0 - MAF[MAF>0.5]
	#y <- y[MAF>0.2,]
	#MAC, made this a lot less
	y <-y[MAF>0.05,]
	# now get the fsts
	fsts <- y$Fst;
	num.snps <- length(fsts)

	pq <- pob[[C]]$Perm$Fst  # permutation quantile matrix.  Just nice to have a short variable

	if(C=="Z") {
		max.y <- 1 #z.chrom.y.max
	} else {
		max.y <- autosome.y.max #oddly, this is NA for me MAC
	}
	new.vals <- plot.chromo.sparkline (
							chrom=C,
							chrom.len=zf.cl$NumBases[ zf.cl$Chromo==C ],
							pos=y$GenCoord,
							fsts=fsts,
							chrom.ht=chrom.height,
							x.span=max.len,
							start.line=start.line,
							x.start=x.start,
							kx=pob[[C]]$Obs$Fst$z,
							ky.list=list(pob[[C]]$Obs$Fst$yhat, pq[, "99.9%"], pq[, "99.99%"], pq[, "100%"], pq[, "0%"]), # list of kernel density y vectors. [[1]] is observed, [[2]] is 1/1000, [[3]] is 1/10000, [[4]] is max of the upper significance stuff 1 in 25K.
#ky.list has NA
							#plot.genes=T,
							max.y=max.y,
							#StartEnds=IslandStartAndEnds[[C]]
#MAC removing references to islands
							)
	start.line <- new.vals[[1]]
	x.start <- new.vals[[2]] + .07 * max.len  # the .1 here dictates how much space we have after the end of the chromo before the start of the next one
	print(new.vals)
}

# now, at the very end we have to draw a scale for megabases
scl <- 131e6 # scale left end point
sch <- 8.34  # scale height
scr <- scl+20e6  # This is how long it will be
scth <- .1 # height of the scale ticks

segments(scl, sch, scr, sch)
segments(seq(scl,scr,1e06), sch, seq(scl,scr,1e06), sch-scth)  # minor ticks at 1 Mb
segments(seq(scl,scr,5e06), sch, seq(scl,scr,5e06), sch-2.0*scth)  # major ticks at 5 Mb
text(seq(scl,scr,5e06), sch-2.0*scth, labels=c(0,5,10,15,20), pos=1, cex=0.75)  # axis tick labels
text((scl+scr)/2, sch+0.5*scth, labels="Horizontal Scale (Mb)", pos=3, cex=0.85)

#dev.off()
dev.copy2pdf(file="MultiChromFig.pdf")
