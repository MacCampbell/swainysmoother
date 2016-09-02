#! /bin/bash
# This script will feed variables to all the other scripts. It needs the following to be specified:
# Working directory, popsize 1, popsize 2, nreps, nprocessors
#usage ./driveSwainySmoother.sh /home/mac/swainysmoother/exampleData/ 12 14 10 5
# for example data ./driveSwainySmoother.sh "/home/mac/swainysmoother/exampleData" 2 2 10 5
# this will make a /home/mac/swainysmoother/exampleData/out directory with output in it

dir=$1; # /home/mac/swainysmoother/exampleData
pop1=$2; # 12
pop2=$3; # 14
nreps=$4; # 10
nprocessors=$5; # 5
#making "out" dir
mkdir $dir/out

# Step 1, run swainyPrepper2.R
Rscript --vanilla ../R/swainyPrepper2.R $dir $pop1 $pop2

# Step 2, run swainyEntryMaker.R
Rscript --vanilla ../R/swainyEntryMaker.R $dir

# Step 3 set up swthFunc.R and perm_parallel.R
for C in "omy01" "omy02" "omy03" "omy04" "omy05" "omy06" "omy07" "omy08" "omy09" "omy10" "omy11" "omy12" "omy13" "omy14" "omy15" "omy16"  "omy17" "omy18" "omy19" "omy20" "omy21" "omy22" "omy23" "omy24" "omy25" "omy26" "omy27" "omy28" "omy29"; do
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


