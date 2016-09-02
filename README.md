# swainysmoother
A method for detecting high Fst genomic regions adapted from Ruegg et al. (2014), doi:10.1111/mec.12842.

## Fst Outliers
This method is to identify genomic regions that are elevated in Fst from RADseq data between a pair of populations. Coverage and SNP density are not uniform in experiments of this type. The method accounts for loci length and SNP density to clearly identify genomic regions in which Fst is elevated. See the original Ruegg et al. (2014) paper for details. Original code by Eric Anderson is available here: http://datadryad.org/resource/doi:10.5061/dryad.73gj4

## Organisms
Originally designed from Swainson's Thrush, I have adapted the code for a chromosome level assembly of rainbow trout. 

## Steps
To use this code yourself, you'll need to make a few changes. Primarily, the number of chromosomes and names need to be changed in numerous locations. The sex chromosomes need some special consideration. Otherwise these are the basic steps.

1. Map your RADseq reads to your reference genome. 
2. Convert these files to .dat files with "processPileup.sh." Take a look at this file, the awk command can be used generally but you'll have to change the samtools command for your specific dataset.
3. Put all these .dat files in a single directory. The .dat files should be in order pop1Sample1.dat pop1Sample2.dat pop2Sample1.dat pop2Sample2.dat as the R script simply reads them in order.
4. from ~/swainysmoother/bashScripts/ you can kick off the analysis by indicating the directory containing your .dat files, the population sizes, the number of reps to run, and the number of processors to use (see driveSwainySmoother.sh for my detail):
~/swainysmoother/bashScripts$ ./driveSwainySmoother.sh "/home/mac/freshwater/combinedBigCreek" 12 14 10000 8
5. Inside the specificied directory, a new subdirectory "out" will have been created containing the necessary files for multiChromoFigure.R and retrieveHighFstSnps.R




