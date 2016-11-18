#! /bin/bash 

#for f in *.sorted.bam;
#	do
#		samtools mpileup -f /home/mac/data/OmyGenomeV06/omyV6Chr.fasta  $f > "`basename $f .sorted.bam`.pileup"
#	done

for f in *.pileup;
	do
		awk -v thresh=5 '$4 >= thresh {if(on==1) {sum+=$4; n++; next;}  start = $2; on = 1; sum=$4; n=1;} $4<thresh {if(on==0) next; on=0; print $1, start, $2, sum/n}' $f > "`basename $f .pileup`.5x.dat"
	done


