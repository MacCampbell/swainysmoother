# get a list of chromosomes in the order that we want them in
#this won't overide an existing .rda file!
args = commandArgs(trailingOnly=TRUE)

#CHR <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","LG2","LGE22")
#CHR <- c("omy01", "omy02", "omy03", "omy04")
CHR <- c("omy01", "omy02", "omy03", "omy04", "omy05", "omy06", "omy07", "omy08", "omy09", "omy10", "omy11", "omy12", "omy13", "omy14", "omy15", "omy16",  "omy17", "omy18", "omy19", "omy20", "omy21", "omy22", "omy23", "omy24", "omy25", "omy26", "omy27", "omy28", "omy29")

REPS <- as.numeric(args[1])  # we need this to form the file names

#setwd("/home/mac/swth/trout/")
#setwd("/home/mac/freshwater/combinedBigCreek/swainy")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis/swainy")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel/swainy")
setwd(args[2])
#load("/Users/tanya/Desktop/Science Fair/swth_snps.rda")
#load("/home/mac/freshwater/combinedBigCreek/data.rda")
#load("/home/mac/freshwater/combinedBigCreekSanLuis/data.rda")
#load("/home/mac/freshwater/combinedBigCreekGabriel/data.rda")
load("./out/data.rda")
for(i in CHR) {
	load(file.path("./out",paste("chr_", i, "_quants_", REPS, ".rda", sep="")))
}


# now, let's also get the observed values and the x values:
source("~/swainysmoother/R/swthFunc.R")


autosomes <- swth.SNPs[swth.SNPs$Chromosome!="Z",]
autosomes <- droplevels(autosomes)
allsnps <- DNAstats(autosomes)
#MAC there are fifteen nans in the trout data...
allsnps[is.na(allsnps)]<-0
allstuff <- smoothval(allsnps, keep.z=TRUE)

# now allstuff is a list of autosomes named correctly and in order.
# so now we can merge all that into a big list.  Maybe we could make data frames.
# Actually, it seems cumbersome to squash them all into data frames.  Let's just keep
# things as a big list.  Each chromosome will have components Obs and Perm which stores
# the results
PermAndObs <-	lapply(CHR,
	function(C) {
		list(Obs=allstuff[[C]], Perm=get(paste("chr_", C, "_quants_", REPS, sep="")))
	}
)
names(PermAndObs) <- CHR

assign(paste("PermAndObs", REPS, sep="_"), PermAndObs)


save(
	list=c(paste("PermAndObs", REPS, sep="_")),
	file=paste("./out/PermAndObs_", REPS, ".rda", sep=""),
	compress="xz"
)
