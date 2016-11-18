#! /bin/bash
# bashScripts$ ./driveSwainySmootherWGS.sh "/home/user/directory" 2 2 25000 10

dir=$1; # /home/user/directory
pop1=$2; # 2
pop2=$3; # 2
nreps=$4; # 25000
nprocessors=$5; # 10
#making "out" dir
mkdir $dir/out

# Step 1, run swainyPrepperWGS.R
Rscript --vanilla ../R/swainyPrepperWGS.R $dir $pop1 $pop2

# Step 2, run swainyEntryMaker.R
Rscript --vanilla ../R/swainyEntryMaker.R $dir

# Step 3 set up swthFunc.R and perm_parallel.R
for C in "omy01" "omy02" "omy03" "omy04" "omy05" "omy06" "omy07" "omy08" "omy09" "omy10" "omy11" "omy12" "omy13" "omy14" "omy15" "omy16"  "omy17" "omy18" "omy19" "omy20" "omy21" "omy22" "omy23" "omy24" "omy25" "omy26" "omy27" "omy28" "omy29"; do
#for C in "omy05"; do

	(
		echo "REPS <- $nreps"
		echo "C<-\"$C\""
		echo "source(\"~/swainysmoother/R/swthFunc.R\")"
		echo "load(\"$dir/out/data.rda\")"
		cat ~/swainysmoother/R/perm_parallel.R
	) > $dir/out/rscript_$C.R
done

# Step 4, set up and run parallel
for i in $dir/out/rscript_*.R; do echo "R CMD BATCH --vanilla $i;"; done > $dir/out/comms4parallel.txt
cd $dir/out/
parallel -j $nprocessors < $dir/out/comms4parallel.txt

# Step 5, compile parallel and compute observed
cd ~/swainysmoother/bashScripts
Rscript --vanilla ../R/compile_parallel_and_compute_observed.R $nreps $dir

# Step 6, basic popgen
Rscript --vanilla ../R/BasicPopGenStatistics.R $dir

# Step 7, make figures (not elegantly done)
Rscript --vanilla ../R/multiChromoFigure.R $dir
Rscript --vanilla ../R/multiChromoFigure-11-20.R $dir
Rscript --vanilla ../R/multiChromoFigure-21-29.R $dir
Rscript --vanilla ../R/multiChromoFigureOmy05.R $dir

# Step 8, retrieve high Fst SNP locations
Rscript --vanilla ../R/retrieveHighFstSnps.R $dir
