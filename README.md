# swainysmoother
A method for detecting high Fst genomic regions adapted from Ruegg et al. (2014), 10.1111/mec.12842. 
[I'm an inline-style link](https://www.google.com)


Applied in [10.1093/jhered/esab049.] (https://academic.oup.com/jhered/advance-article/doi/10.1093/jhered/esab049/6356224?searchresult=1#:~:text=https%3A//doi.org/10.1093/jhered/esab049)



## Fst Outliers
This method is to identify genomic regions that are elevated in Fst from RADseq data between a pair of populations. Coverage and SNP density are not uniform in experiments of this type. The method accounts for loci length and SNP density to clearly identify genomic regions in which Fst is elevated. See the original Ruegg et al. (2014) paper for details. Original code by Eric Anderson is available here: http://datadryad.org/resource/doi:10.5061/dryad.73gj4

## Organisms
Originally designed from Swainson's Thrush, I have adapted the code for a chromosome level assembly of rainbow trout. 

## Steps
To use this code yourself, you'll need to make a few changes. Primarily, the number of chromosomes and names need to be changed in numerous locations. The sex chromosomes need some special consideration. Otherwise these are the basic steps.

1. Map your RADseq reads to your reference genome creating sam and sorted bam files.
2. Call genotypes, and make an "out.geno" file. I used angsd like so:

angsd -bam combinedPop -out genotypes -doSaf 1 -doPost 2 -GL 2 -doMaf 2 -doMajorMinor 1 -minMaf 0.05 -minInd 22 -minMapQ 30 -minQ 30 -SNP_pval 2e-3 -anc ../omyV6Chr.fasta  -postCutoff 0.95 -doGeno 2
gunzip -c genotypes.geno.gz | cut -f 1-28 > out.geno


3. Convert alignment files to .dat files with "processPileup.sh." Take a look at this file, the awk command can be used generally but you'll have to change the samtools command for your specific dataset.
4. Put all these .dat files in a single directory. The .dat files should be in order pop1Sample1.dat pop1Sample2.dat pop2Sample1.dat pop2Sample2.dat as the R script simply reads them in order.
5. A file of chromosome lengths (RT.Chrom.Lengths) needs to be in the same length directory as the .dat and out.geno files. One way to get the lengths of your chromosomes is: while read -r line; do echo ${#line}; done < omyV6Chr.fasta >> chromLengths.txt. A properly formatted file is in exampleData.
6. From ~/swainysmoother/bashScripts/ you can kick off the analysis by indicating the directory containing your .dat files, the population sizes, the number of reps to run, and the number of processors to use (see driveSwainySmoother.sh for my detail):
~/swainysmoother/bashScripts$ ./driveSwainySmoother.sh "/home/mac/freshwater/combinedBigCreek" 12 14 10000 8
7. Inside the specificied directory, a new subdirectory "out" will have been created containing the necessary files for multiChromoFigure.R and retrieveHighFstSnps.R

##Whole Genome Resequencing Data
Whole genome resequencing data is not continuous in the amount of coverage over a genome. To work with 10X coverage whole genome resequencing data, I have altered some of the swainysmoother scripts. The original processPileup.sh filters for 10X, and with the data I worked with this did not produce many overlaps of any length between individuals. Therefore, I made a processPileup5x.sh file (see exampleWGS directory). Furthermore, the swainysmoother "swainyPrepper" script will spend a long time looking for overlap with very small fragments of continuous coverage, so I put together a simple bash/baseR script pair to filter .dat files (exampleWGS directory, subDat.sh, sub.R). I examined the length of segments of continuous coverage with the "evaluateDat" scripts in exampleWGS. These steps are going to be fairly unique to each case examined. I also created a "driveSwainySmootherWGS.sh" script to run "swainyPrepperWGS.R" and produce results with the 5x filtered data.



