REPS <- 10
C<-"omy14"
source("~/swainysmoother/R/swthFunc.R")
load("/home/mac/swainysmoother/exampleData/out/data.rda")
# This is a skeletal file.  Some values must be set before
# running this. Namely, C, must be set to the name of a chromosome
#(i.e. something in
# 1 1A 1B 2 3 4 4A 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 LG2 LGE22)
# and REPS must be set to the number of permutations you want to do


# read in necessary functions
#source("/Users/Mac/rstudio/swainysmoother/trout/swthFunc.R")
#source("/home/mac/swainysmoother/R/swthFunc.R")

# read in the data

#load("/Users/eriq/Documents/work/teaching/indep_study/jazz_pouls_2013/data/swth_snps_clean.rda")
#load("/Users/tanya/Desktop/Science Fair/swth_snps.rda")
#load("/Users/Mac/rstudio/swainysmoother/trout/data.rda")
#load("/home/mac/freshwater/combinedBigCreek/data.rda")
#load("/home/mac/freshwater/combinedBigCreekSanLuis/data.rda")
#load("/home/mac/freshwater/combinedBigCreekGabriel/data.rda")
#load("/home/mac/freshwater/combinedBigCreekMal/data.rda")

# drop the Z chromosome from the data set
autosomes <- swth.SNPs[swth.SNPs$Chromosome!="Z",]
autosomes <- droplevels(autosomes)
allsnps <- DNAstats(autosomes)
#again, some NA's show up
allsnps[is.na(allsnps)]<-0

# split all the data into chromosomes
in.chroms <- split.and.aggregate(allsnps)

# grab as many SNPs at random (without replacement) from the entire set of autosomes
# as needed to fill up the SNPs in chromosome C with they randomly drawn values
ps <- lapply(1:REPS, function(x) {smoothval.single.chrom(permSNPsSingleChrom(in.chroms$chr[[C]],allsnps))})


# now we want to get the quantiles off those.  We need to do it for each of the
# different statistics we are tracking:
stats <- list(df="df", dxy="dxy", Pi.1="Pi.1", Pi.2="Pi.2", Fst="Fst")

# here is the variable name we would like to assign the result to
var.name <- paste("chr", C, "quants", REPS, sep="_")

# then we use "assign" because we want to assign it to a variable name
# that we made up using the variables C and REPS
#MAC some errors here na.rm=TRUE
assign(var.name,
	value=lapply(stats, function(st) {
	t(apply(sapply(ps, function(x) x[[st]]$yhat), 1, quantile, na.rm=TRUE, probs=c(0, 1e-4, 1e-3, 1e-2, .05, .25, .5, .75, .95, 1-1e-2, 1-1e-3, 1-1e-4,1 )))
	})
)


# then we save that variable to an rda file.  Note that the symbol is the value of var.name:
save(list=c(paste(var.name)), file=paste(var.name,".rda", sep=""), compress="xz")
#save(list=c(paste(var.name)), file=paste($dir,"/out/",var.name,".rda", sep=""), compress="xz")
