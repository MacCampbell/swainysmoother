# swainysmoother
A method for detecting high Fst genomic regions adapted from Ruegg et al. (2014), doi:10.1111/mec.12842.

## Fst Outliers
This method is to identify genomic regions that are elevated in Fst from RADseq data between a pair of populations. Coverage and SNP density are not uniform in experiments of this type. The method accounts for loci length and SNP density to clearly identify genomic regions in which Fst is elevated. See the original Ruegg et al. (2014) paper for details. Original code by Eric Anderson is available here: http://datadryad.org/resource/doi:10.5061/dryad.73gj4

## Organisms
Originally designed from Swainson's Thrush, I have adapted the code for a chromosome level assembly of rainbow trout. 
