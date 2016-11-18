#! /bin/bash

for f in ./originalDat/*.dat;

	do echo "`basename $f .5x.dat`.sub.dat";
	./sub.R $f "sub1"
	cut -f 2-5 -d ' ' "sub1" |  perl -pe 's/"//g' > sub2
	#remove the first line
	echo "$(tail -n +2 sub2)" > "`basename $f .5x.dat`.sub.dat"
	"yes" | rm sub1
	"yes" | rm sub2
done;
