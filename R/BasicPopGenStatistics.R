# this is script that computes statistics for each SNP and also
# aggregates those into contig specific measures as needed

# first get the data:
# Set to your own working directory:
#setwd("/Users/eriq/Documents/work/assist/kristenruegg/SWTH_genome_analysis_for_paper")
#setwd("/home/mac/freshwater/combinedBigCreek")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis/")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel")
setwd("/home/mac/freshwater/combinedBigCreekMal")
# load the data we need:
#load("inputs_and_scripts/swth_snps_mapped_to_ZF.rda")
load("data.rda")

#Mac, kept naming the same

rownames(ZF.Chrom.Lengths) <- ZF.Chrom.Lengths$Chromo  # give these some rownames




########### NOW  WE DEFINE A BUNCH OF FUNCTIONS
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
	as.numeric(xor(x$Pop1.f1 + x$Pop2.f2==0, x$Pop2.f1 + x$Pop1.f2==0))
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
	sumagg <- aggregate(cbind(Pi.1, Pi.2, dxy, df) ~ CHROM, x, sum)
	meanagg <- aggregate(cbind(length, GenCoord) ~ CHROM, x, mean)
	x <- merge(sumagg,meanagg)
	x <- x[order(x$GenCoord),]
	x
}

# this function adds calculated variables to each SNP
DNAstats <- function(x) {
 	#MAC length is now length
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
############## END OF FUNCTIONS

# split data into a list of chromosomes and drop levels from factors in each
chr <- split(swth.SNPs, swth.SNPs$Chromosome, drop=TRUE)
lapply(chr, droplevels) -> chr


# get the stats for each SNP in the whole data set\
SNP.stats <- lapply(chr, DNAstats)

# reordering names
chr.name <- gsub("B", ".2", gsub("A", ".1", names(SNP.stats)))
#chr.name[31:33] <- 31:33
#SNP.stats <- SNP.stats[order(as.numeric(chr.name))]
SNP.stats <- SNP.stats[order(chr.name)]
contig.stats <-lapply(SNP.stats, agg.contig)


# now save those variables in an Rda we can load later as needed
save(SNP.stats, contig.stats, file="SNP_and_Contig_Stats.rda", compress="xz")



# this is other remnant stuff that was done for the permutation bits
#autosomes <- swth.SNPs[swth.SNPs$Chromosome!="Z",]
#autosomes <- droplevels(autosomes)

#allsnps <- DNAstats(autosomes)
#chr1000 <- split(allsnps, allsnps$Chromosome, drop=T)
#chr.name <- gsub("B", ".2", gsub("A", ".1", names(chr1000)))
# chr.name[31:33] <- 31:33
#chr1000 <- chr1000[order(as.numeric(chr.name[1:32]))]

#contig <-lapply(chr1000, agg.contig)


