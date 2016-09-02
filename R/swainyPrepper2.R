#This program reads in filtered pilup files from a directory
#On these files, the will be made into a windows and population genetic statistics computed
#starting with *.dat files

#library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
args = commandArgs(trailingOnly=TRUE)

#set your working directory
#setwd("/home/mac/freshwater/combinedBigCreek")
#setwd("/home/mac/freshwater/combinedBigCreekSanLuis")
#setwd("/home/mac/freshwater/combinedBigCreekGabriel")
#setwd("/home/mac/freshwater/combinedBigCreekMal")

setwd(args[1])
#The size of each population needs to be specified
popSize1<-as.numeric(args[2])
popSize2<-as.numeric(args[3])

#popSize1<-12 #bigCreek comparison
#popSize2<-14

#popSize1<-14
#popSize2<-8 #San Luis Rey comparison

#popSize1<-14
#popSize2<-10 #San Gabriel

#popSize1<-14
#popSize2<-18 # Mal

#we need a genotype file, as from ANGSD, be careful because there may be an extra tabular column
genos<-read.table("./out.geno", sep="\t")

#http://stackoverflow.com/questions/4862178/remove-rows-with-nas-in-data-frame
delete.na <- function(DF, n=0) {
  log <- apply(df, 2, is.na)
  logindex <- apply(log, 1, function(x) sum(x) <= n)
  df[logindex, ]
}

fileNames<-dir("./",pattern=".dat")

halfPop<-length(fileNames)/2
popSize<-length(fileNames)



#files<-lapply(fileNames, read.table, header=FALSE)
#d1<-lapply(files, tbl_df)

df<-NULL
myDiffs<-NULL
myLefts<-NULL
myRights<-NULL
i<-1

#create the first dataframe for later joining
leftCoor<-paste("left",i,sep="")
rightCoor<-paste("right",i,sep="")
diff<-paste("diff",i,sep="")
myDiffs[i]<-diff
myLefts[i]<-leftCoor
myRights[i]<-rightCoor
data <- read.table(fileNames[i],header=FALSE)
d <- tbl_df(data)
dd <- d %>% mutate(trunc1 = (V2 %/% 10000) * 10000, trunc2 = (V3 %/% 10000) * 10000, length = V3 - V2)

ddd <- dd %>% arrange(desc(length)) %>% group_by(V1, trunc1, trunc2) %>% summarise(left = min(V2), right = max(V3)) %>% mutate(diff = right - left) %>% filter(diff >150 & diff < 600)
ddd<-ddd %>% rename(chrom=V1)
names(ddd)[4]<-leftCoor
names(ddd)[5]<-rightCoor
names(ddd)[6]<-diff
df<-ddd



for(i in 2:length(fileNames)) {
  #print(i)
  leftCoor<-paste("left",i,sep="")
  #print(leftCoor)
  rightCoor<-paste("right",i,sep="")
  #print(rightCoor)
  diff<-paste("diff",i,sep="")
  myDiffs[i]<-diff
  myLefts[i]<-leftCoor
  myRights[i]<-rightCoor
  data <- read.table(fileNames[i],header=FALSE)
  d <- tbl_df(data)
  dd <- d %>% mutate(trunc1 = (V2 %/% 10000) * 10000, trunc2 = (V3 %/% 10000) * 10000, length = V3 - V2)

  ddd <- dd %>% arrange(desc(length)) %>% group_by(V1, trunc1, trunc2) %>% summarise(left = min(V2), right = max(V3)) %>% mutate(diff = right - left) %>% filter(diff >150 & diff < 600)
  #ddd <- dd %>% arrange(desc(length)) %>% group_by(V1, trunc1, trunc2) %>% rename_("left" = "V2", "right" = "V3") %>% mutate(diff = right - left) %>% filter(diff >150 & diff < 600)

  #These comments sometimes work and sometimes don't, I don't understand why. Line 68 and 75 that is. I think it has to do with plyr and dplyr
  #ddd <- dd %>% arrange(desc(length)) %>% group_by(V1,trunc1,trunc2) %>% mutate(diff = V3 - V2) %>% filter(diff > 100 & diff < 600)
  #ddd <- dd %>% arrange(desc(length)) %>% group_by(V1, trunc1, trunc2) %>% rename_("left" = "V2", "right"="V3") %>% mutate(diff = right - left) %>% filter(diff >150 & diff < 600)

  #ddd<-plyr::rename(ddd, c(V1="chrom", "left"=leftCoor, "right"=rightCoor,"diff"=diff))
  ddd<-ddd %>% rename(chrom=V1)
  names(ddd)[4]<- leftCoor
  names(ddd)[5]<-rightCoor
  names(ddd)[6]<-diff
   #ddd<-plyr::rename(ddd, c(V1="chrom", V2=leftCoor,V3=rightCoor, "diff"=diff))
    assign(paste("table",i,sep=""),ddd)

    #df<-merge(df,ddd,all=TRUE)
    df<-full_join(df,ddd)
}

dft <- tbl_df(df)
#get a unique list of chromosomes
chromos<-unique(dft$chrom)
#extract from dft those entries with data from at least 1/2 popsize for all diff columns
#there are three entries, so 3 x NA threshold, left, right, diff
dfClean<-delete.na(df, (3 * halfPop))
dfC<-tbl_df(dfClean)
#now we need min/maxes for each row for left/right bounds
dfC$leftEnd<-apply(dfC[myLefts],1,min,na.rm=TRUE)
dfC$rightEnd<-apply(dfC[myRights],1,min,na.rm=TRUE)
#and make a length column that has the maximum extent of the left and right bounds
dfCC <- dfC %>% mutate(length=(rightEnd - leftEnd) )


#getting genotype file, finding overlap

genos<-tbl_df(genos)
#edited 8/13/2016

genos<-rename(genos, chrom=V1, site=V2)
#genos<-plyr::rename(genos, c(V1="chrom", V2="site"))

myLoci<-NULL
myLoci$chrom<-dfCC$chrom
myLoci$left<-dfCC$leftEnd
myLoci$right<-dfCC$rightEnd
myLoci$length<-dfCC$length
myLoci<-as.data.frame(myLoci)
myLoci<-tbl_df(myLoci)
#now to filter for overlap (takes a minute)
#produces a chrom site and which contig a site is on
contigs <- genos %>% group_by(chrom, site) %>% summarise(contig_idx = paste(which(myLoci$left < site & myLoci$right > site), collapse = ","))
myLoci$contig_idx <- 1:nrow(myLoci)
#fulls has contigs with which sites are in it
fulls <- full_join(myLoci, contigs %>% mutate(contig_idx = as.integer(contig_idx)))
#fulls allows us to identify if there are any sites in multiple contigs

#filteredContigs<-contigs %>% filter(str_detect(contig_idx, ",")) %>% left_join(., myLoci)
#filteredContigs<-full_join(myLoci, contigs %>% mutate(contig_idx = as.integer(contig_idx)))
#remove NA entries from fulls and remove duplicated sites
fulls2<-fulls %>% group_by(site) %>% filter(row_number() == 1)
#don't want NA's
fulls3<-fulls2 %>% filter(site >= 1)

#Let's join full3 with genos, has redundant chrom column now...
filteredData<-right_join(fulls3, genos, by = c("site","chrom")) %>% filter(contig_idx >= 1)
#genotypes, calculate allele occurences
subseted<-filteredData[,grepl("V",colnames(filteredData))]

#now we need to subset and get freqs for each pop

pop1Vector<-NULL
pop2Vector<-NULL
for(i in 1:popSize1) { pop1Vector[i]<-paste("V",i+2,sep="")}
for(i in 1:popSize2) { pop2Vector[i]<-paste("V",i+2+popSize1,sep="")}
freqs1<-NULL
freqs2<-NULL
subset1<-filteredData[pop1Vector]
subset2<-filteredData[pop2Vector]

freqs1<-NULL
freqs2<-NULL
freqs1$site<-filteredData$site
freqs2$site<-filteredData$site
freqs1$chrom<-filteredData$chrom
freqs2$chrom<-filteredData$chrom

freqs1$zero<-rowSums(subset1==0)
freqs1$one<-rowSums(subset1==1)
freqs1$two<-rowSums(subset1==2)
freqs1$missing<-rowSums(subset1==-1)
freqs1df<-tbl_df(freqs1)
freqs2$zero<-rowSums(subset2==0)
freqs2$one<-rowSums(subset2==1)
freqs2$two<-rowSums(subset2==2)
freqs2$missing<-rowSums(subset2==-1)
freqs2df<-tbl_df(freqs2)

calcFreqs1<-freqs1df %>% mutate(Pop1.f1 = ((2 * zero) + (one)), Pop1.f2=((2 * two)+(one)))
calcFreqs2<-freqs2df %>% mutate(Pop2.f1 = ((2 * zero) + (one)), Pop2.f2=((2 * two)+(one)))
subCalcFreqs1<-select(calcFreqs1, chrom, site, Pop1.f1, Pop1.f2)
subCalcFreqs2<-select(calcFreqs2, chrom, site, Pop2.f1, Pop2.f2)
totCalcFreqs<-left_join(subCalcFreqs1,subCalcFreqs2, by = c("chrom","site"))
allData<-full_join(filteredData, totCalcFreqs, by=c("chrom","site"))
#get what we really want-chrom, left, right, length, contig_idx, site, f1, f2
final<-tbl_df(select(allData, chrom, left, right, length, contig_idx, site, Pop1.f1, Pop1.f2, Pop2.f1, Pop2.f2))
filteredFinal<-final %>% na.omit
#or filter(!is.na(left))
#desired output for one pop
write.table(filteredFinal,file="./out/convertedData.txt")
