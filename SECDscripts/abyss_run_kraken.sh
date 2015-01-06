#!/bin/sh

# abyss_run.sh
# Purpose: to run abyss-pe on R1 and R2 files
# contained in Zips directory

# Set directory to copy assembly output
#abyssFile="/Volumes/Data_HD/Brucella/ceti/abyss_files"

# Change working directory to directory containing Read Files
mkdir ./Zips
mv *fastq* ./Zips

cd ./Zips
if [ -f *R2* ]; then
    echo "R2 paired end read file present"
    forReads=`ls | grep _R1`
    echo "Forward Reads:  $forReads"
    revReads=`ls | grep _R2`
    echo "Reverse Reads:  $revReads"
else
    echo "Just a single read present"
    forReads=`ls | grep *fastq`
    echo "Forward Reads:  $forReads"
fi


# Place R1 Reads into variable
#forReads=`ls | grep _R1`
#echo "Forward Reads:  $forReads"
# Place R2 Reads into variable
#revReads=`ls | grep _R2`
#echo "Reverse Reads:  $revReads"

# Get name of isolate to use for naming output files
n=`echo $forReads | sed 's/_.*//' | sed 's/\..*//'`
echo "***Isolate naming convention:  $n"

# De Novo assembly using ABySS
if [ -f *R2* ]; then
    abyss-pe name=${n}_abyss k=64 in="$forReads $revReads"
else
    abyss-pe name=${n}_abyss k=64 in="$forReads"
fi

# Copy contig file to set directory
#cp ${n}_abyss-3.fa $abyssFile

# Clean-up output
# Move reads to their own directory
#rm coverage.hist
mkdir ../${n}_abyss
mv ${n}_abyss* ../${n}_abyss
mv ${n}_abyss* ../${n}_abyss
mv coverage.hist ../${n}_abyss

echo "*** Done ***"

#
#  Created by Stuber, Tod P - APHIS on 09/20/13.
#
