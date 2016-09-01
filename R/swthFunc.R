#need to set path to wksmooth.so
#check Mac comments to edit number of chromosomes

#### Functions to Compute SNP and contig attributes  ##########
# this function computes SNP's contribution to pi in contig
pifunc <- function(ref, alt, length) {
	(2/(length*(ref+alt)*((ref+alt)-1))) * ref * alt
	 }


# this function computes SNP's contribution
# to pairwise nucleotide diveristy in interspecific comparisons
dxycomp <- function(x) {
	((x$Pop1.f1 * x$Pop2.f2) + (x$Pop2.f1 * x$Pop1.f2)) /
	((x$Pop1.f1 + x$Pop1.f2) * (x$Pop2.f1 + x$Pop2.f2)*x$length)
}


# this function determines if the SNP is a fixed difference or not
dfcomp <- function(x) {
	as.numeric(xor((x$Pop1.f1 + x$Pop2.f2==0), (x$Pop2.f1 + x$Pop1.f2==0)))
}

# this function computes Fst, based on Pi.1, Pi.2, Pi.All
fstcomp <- function(x) {
	n1 <- x$Pop1.f1 + x$Pop1.f2
	n2 <- x$Pop2.f1 + x$Pop2.f2
 	1 - (  (x$Pi.1*n1*(n1-1)/2  +  x$Pi.2*n2*(n2-1)/2) /
 				(x$Pi.All*((n1*(n1-1)/2)  +  (n2*(n2-1)/2)))  )
 }

# function creates dataframe with 1 entry for each contig in a chromosome
agg.contig <- function(x) {
	#MAC we now have contig_idx, renamed this to CHROM
  sumagg <- aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ CHROM, x, sum)
	meanagg <- aggregate(cbind(length, GenCoord) ~ CHROM, x, mean)
	x <- merge(sumagg,meanagg)
	x <- x[order(x$GenCoord),]
	x
}

# this function adds calculated variables
DNAstats <- function(x) {
#MAC length is now numeric so...
	x$length
	x$Pi.1 <- pifunc(x$Pop1.f1, x$Pop1.f2, x$length)
	x$Pi.2 <- pifunc(x$Pop2.f1, x$Pop2.f2, x$length)
	x$Pi.All <- pifunc(x$Pop1.f1+x$Pop2.f1, x$Pop1.f2+x$Pop2.f2, x$length)
	x$dxy <- dxycomp(x)
	x$df <- dfcomp(x)
	x$Fst <- fstcomp(x)
	x <- x[order(x$GenCoord),]
	x
}



##### FUNCTIONS TO PERMUTE SNPS HOLDING LOCATIONS STEADY
# here is a function to permute the SNP values amongst the SNP locations
permSNPs <- function(x) {
	cbind(x[,c(1:5,12)], x[sample(nrow(x),nrow(x)),-c(1:5,12)])
}


# given a data frame for  single chromosome, insert permuted values from
# the entire genome (in a data frame x) into it
permSNPsSingleChrom <- function(c,x) {
	cbind(c[,c(1:5,12)], x[sample(nrow(x),nrow(c)),-c(1:5,12)])
}




# WKSMOOTH STUFF
# this is a function to do a weighted sliding window smooth.
# It calls a C function that must already be
# in a dynamically loaded shared object library

#dyn.load("/Users/tanya/Desktop/Science Fair/wksmooth.so")
#dyn.load("/Users/Mac/rstudio/swainysmoother/src/wksmooth.so")
dyn.load("/home/mac/swth/src/wksmooth.so")

wksmooth <- function(
	G,  # x-positions of the data
	y,  # measured y-values at each G
	L,  # the weights that apply to each y
	wind.width=(max(G,na.rm=T) - min(G, na.rm=T))/200,  # sliding window width.  By default 1/200th of the range of data
	num.pts=length(G)  # number of equally space points within the range of G at which to compute smoothed values
	)
{
	# check for errors and things
	if(length(y) != length(G)) stop("Error! y and G not of same length!!")
	if(length(L) != length(G)) stop("Error! L and G not of same length!!")
	# we could check for other stupid errors like windows that are obviously too wide, but we won't!

	# make sure things are in sorted order by G
	ord<-order(G)
	G<-G[ord]
	y<-y[ord]
	L<-L[ord]

	z=seq(from=G[1], to=G[length(G)], length.out=num.pts)  # set of points at which to compute smoothed values

	# now, make the call
	#MAC
	#	boing <- .C("wksmooth",

	boing <- .C("wksmooth",
		as.double(G),
		as.double(y),
		as.double(L),
		as.double(z),
		as.integer(length(G)),
		as.integer(length(z)),
		as.double(wind.width),
		result=double(length=length(z))
		)

		list(z=z, yhat=boing$result, G=G, y=y)

}
#Alter here to identify actual coordinates?
wk.wrap <- function(dframe, a, keep.z=F) {
	 wksm <- wksmooth(dframe$GenCoord, dframe[[a]],
	 				  dframe$length, wind.width=1e06)
		if(keep.z==F) {
			list(yhat=wksm$yhat)
		} else  {
			list(z=wksm$z, yhat=wksm$yhat)
		}
}







islands <- function(a,b) {
	n <- 0	# value of divergence index
	L <- length(a)	# length of each vector (in this case, 3)
	z <- outer(a,b, FUN = function(x,y) {abs(x-y)})
	while( {zz <- which.min(z); length(zz)>0} ) {
		c <- (zz-1)%/%L + 1 # find column of min
		r <- (zz-1)%%L + 1  # find row of min
		m <- z[zz]	# distance to nearest peak
		n <- (n + m)
		z[r,] <- NA
		z[,c] <- NA

	}
	n
}




##### HIGH-LEVEL  FUNCTIONS TO DO RELATIVELY COMPLETE ANALYSES

split.and.aggregate <- function(x) {
	# Split into list with components for each chromosome
	cc <- split(x, x$Chromosome, drop=T)
	chr.name <- gsub("B", ".2", gsub("A", ".1", names(cc)))
 	#MAC commented out
	#chr.name[31:33] <- 31:33
#MAC only have 29 trout chromos
 	#cc <- cc[order(as.numeric(chr.name[1:32]))]
 	#cc <- cc[order(as.numeric(chr.name[1:29]))]
	cc <-cc[order(chr.name[1:29])]
	# Aggregate across contigs
	contig <-lapply(cc, agg.contig)
	list(chr=cc, contigs=contig)
}


# This function takes all the SNPs with values computed for them (may or may not
# have been permuted) and returns the values of wksmooth for all the variables and
# for all the chromosomes
smoothval <- function(x, keep.z=F) {
	# Split into list with components for each chromosome
	cc <- split(x, x$Chromosome, drop=T)
	chr.name <- gsub("B", ".2", gsub("A", ".1", names(cc)))
 	chr.name[31:33] <- 31:33
#MAC only have 29 trout chromos
 	#cc <- cc[order(as.numeric(chr.name[1:29]))]
	cc <-cc[order(chr.name[1:29])]
	# Aggregate across contigs
	contig <-lapply(cc, agg.contig)

	smooths <- lapply(contig, function(chr) {
			list(df=wk.wrap(chr, "df", keep.z),
				dxy=wk.wrap(chr, "dxy", keep.z),
				Pi.1=wk.wrap(chr, "Pi.1", keep.z),
				Pi.2=wk.wrap(chr, "Pi.2", keep.z))
		})
		for(i in names(cc)) {
			smooths[[i]]$Fst=wk.wrap(cc[[i]], "Fst", keep.z)
		}

		smooths
}

# pass this a data frame of snps for a single chromosome
smoothval.single.chrom <- function(x, keep.z=F) {
	chr <- agg.contig(x)

	list(df=wk.wrap(chr, "df", keep.z),
		dxy=wk.wrap(chr, "dxy", keep.z),
		Pi.1=wk.wrap(chr, "Pi.1", keep.z),
		Pi.2=wk.wrap(chr, "Pi.2", keep.z),
		Fst=wk.wrap(x, "Fst", keep.z))

}
