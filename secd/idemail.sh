#!/bin/sh

# Script will send out a summary email of stats output from idvirus.sh

if [ -z $1 ]; then 
	echo "

echo "Start Time: `date`" > /scratch/report/dailyTime
starttime=`date +%s`

for i in *.fastq*; do 
	n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
	echo "n is : $n"
	mkdir -p $n
	mv $i $n/
 done

echo 'currentdir=`pwd`
for f in *; do 
	cd $currentdir
	echo $f
	cd ./$f 
	idvirus.sh $1 

done'

echo "" >> /scratch/report/dailyTime
echo "End Time: `date`" >> /scratch/report/dailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "Run time: $runtime seconds" >> /scratch/report/dailyTime


# created 2015-08-19 stuber
