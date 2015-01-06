#!/bin/sh

################################################################################################################################
#|||||||||||||||||||||||||||||||||||||||||| FUNCTION:  PARSE BLAST XML FROM VELVET CONTIGS |||||||||||||||||||||||||||||||||||||
################################################################################################################################

NR_CPUS=50

function parseXML () {

abyss-parseXML-Blast5.py ${n}_contig-blastResults-5.xml | sed 's/\[.*\]//g' | sort -k1,1 -k2,2 | awk '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' | awk '{k=$1; a[k]++; b[k]=$0}; END{for (k in a) print a[k], b[k]}' | sort -rnk1,1 | awk 'BEGIN{print "n", "acc", "length", "score", "e-value", "ID"} {print $0}' | column -t > BLAST-summary-${n}-contigIDs.txt

abyss-parseXML-Blast5.py ${n}_contig-blastResults-5.xml | sed 's/\[.*\]//g' | sort -k1,1 -k2,2 | column -t > BLAST-all-${n}-contigIDs.txt

}

#date > log.txt
#logfile="log.txt"

m=`basename $1`; n=`echo $m | sed 's/_.*//' | sed 's/\..*//'`
echo "***Sample naming convention:  ${n}"

cp $1 ${n}-contigsoriginal.fa

# If contigs going in are not all on the same line then newlines are removed to put fasta on a single line, but the name header name is changed.
#awk '{if ($0 ~ />/ ) {print ">"} else print $0}' $1 | tr -d "\n" | sed -e 's:>:\n>\n:g' | awk '{if ($0 ~  />/ ) {print ">contigs-" x++} else print $0}' | grep -v "^$" > ${n}-contigsoriginal.fa

# Remove everything behond the space.  Meant for ">name otherinfo" to just ">name"
sed 's/ .*$//' ${n}-contigsoriginal.fa > ${n}-contigs.fa

# Cut down contigs and this also should put all fasta contigs onto a single line.
sed -e 's/$/=/g' ${n}-contigs.fa | tr -d "\n" | tr "=" "\n" | grep -v "^$" | awk '{ if ($0 ~ /^>/ ) { print $0 } else if (length($O) > 600 ) {print substr($0, 50, 500) } else {print $0}}' > ${n}-contigs3.fa
echo "##### Blasting contig file #####"
echo "Start Time: `date`"

#head doesn't work if sequence is not in single line
#head -20000 contigs3.fa > contigs2.fa
cat ${n}-contigs3.fa > ${n}-contigs2.fa

blastContigs=`grep -c ">" ${n}-contigs2.fa`
echo "$blastContigs BLASTed contigs"

blastn -query ${n}-contigs2.fa -db /usr/local/bin/ncbi-blast-2.2.29+/db/nt -task blastn -num_threads ${NR_CPUS} -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1

parseXML
wait

rm ${n}-contigsoriginal.fa
#rm ${n}-contigs2.fa
mv ${n}-contigs2.fa BLAST-${n}-INFILE.fasta
rm ${n}-contigs3.fa
rm ${n}-contigs.fa
#rm ${n}_contig-blastResults-5.xml

# From xml get individual read headers with identification
egrep "<Iteration_query-def>|<Hit_accession>|<Hit_def>" ${n}_contig-blastResults-5.xml | sed 's/<Iteration_query-def>//g' | sed 's:</Iteration_query-def>:TABPLACEMENT:g' | sed 's/<Hit_def>//g' | sed 's:</Hit_def>:TABPLACEMENT:g' | sed 's:<Hit_accession>::g' | sed 's:</Hit_accession>:NEWLINEPLACEMENT:g' | tr -d "\n" | sed -e 's:NEWLINEPLACEMENT:\n:g' | sed -e 's:TABPLACEMENT:\t:g' > BLAST-headersall-${n}-contigIDs.txt

## Possibly look at BLAST summary and grep Accessions to isolate all reads.

## Get reads from header file
#grep -F -A3 -h -f (*Header File*) (*Raw Fastq*) | grep -v '^--$' >> IsolatedReads.fastq

## Relabel abyss headers
#awk '{if ($0 ~  />/ ) {print ">abyss-" x++} else print $0}' Abyss-8-edit1.fasta > Abyss-8-edit2.fasta
