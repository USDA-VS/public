#!/bin/sh

# Set varibles
test=/bioinfo11/MKillian/Analysis/script_dependents/ai/testflu/*fasta
root=`pwd`
sampleName=`ls *.fastq* | head -1 | sed 's/_.*//' | sed 's/\..*//'`

function findbest () {

ref=`ls | grep .fasta`
echo "Reference Input:  $ref"
refname=${ref%.fasta}

bwa index $ref
samtools faidx $ref
java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${refname}.dict
    if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allow
if [[ $sampleType == "paired" ]]; then
    bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"$refname""\t"PL:ILLUMINA"\t"PU:"$refname"_RG1_UNIT1"\t"LB:"$refname"_LIB1"\t"SM:"$refname" $ref $forReads $revReads > ${refname}.sam
else
    bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"$refname""\t"PL:ILLUMINA"\t"PU:"$refname"_RG1_UNIT1"\t"LB:"$refname"_LIB1"\t"SM:"$refname" $ref $forReads > ${refname}.sam
fi

samtools view -bh -F4 -T $ref ${refname}.sam > ${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${refname}.raw.bam ${refname}.sorted
echo "****Indexing Bam"
samtools index ${refname}.sorted.bam

#Number of nucleotides in reference with coverage
echo "*** Bamtools is getting coverage..."
bamtools coverage -in ${refname}.sorted.bam | awk -v x=${refname} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${refname}-coveragefile

#Mean depth of coverage
meancov=`awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' ${refname}-coveragefile`

#Length of reference
countNTs=`awk 'END{print NR}' ${refname}-coveragefile`

#count positions with coverage
covCount=`awk '{ if ($3 != 0) count++ } END { print count }' ${refname}-coveragefile`
echo "covCount $covCount"

declare -i x=${covCount}
declare -i y=${countNTs}

#Percent of reference with coverage
perc=`awk -v x=$x -v y=$y 'BEGIN { print(x/y)*100}'`
echo "perc: $perc"

printf "%-20s %11.2f%% %'10dX\n" ${refname} $perc $meancov >> ${root}/${s}/${sampleName}.findbest
}

#Establish Read Files
if [ -f *R2* ]; then
    echo "R2 paired end read file present"
    export sampleType="paired"
    forFile=`ls | grep _R1`
    forReads="$root/$forFile"
    echo "Forward Reads:  $forReads"
    revFile=`ls | grep _R2`
    revReads="$root/$revFile"
    echo "Reverse Reads:  $revReads"
else
    echo "Just a single read present"
    export sampleType="single"
    forFile=`ls | grep fastq`
    forReads=$root/$forFile
    echo "Forward Reads:  $forReads"
fi

#if [ "$bflag" ]; then
    echo ""
    echo " *** B FLAG ON, BUG FINDING MODE, SINGLE SAMPE PROCESSING *** "
    echo ""
    cp $test ./
	cd ${root}
	mkdir segment{1..8}
	mv segment1*fasta ./segment1/
        mv segment2*fasta ./segment2/
	mv segment3*fasta ./segment3/
	mv segment4*fasta ./segment4/
	mv segment5*fasta ./segment5/
	mv segment6*fasta ./segment6/
	mv segment7*fasta ./segment7/
	mv segment8*fasta ./segment8/
	for s in segment*; do
		cd $root
		cd $s
		echo "s: $s"
		pwd
		for i in *fasta; do
			(cd ${root}/${s}
			mkdir ${i%.fasta}
        		mv ${i} ${i%.fasta}
			ln ${root}/*fastq* ${i%.fasta}
        		echo "working on $sampleName $s $i"
       			cd ${i%.fasta}; findbest) &
	        let count+=1
                [[ $((count%55)) -eq 0 ]] && wait
    		done
	wait
	cd ${root}/$s
	best=`sort -rk2,3 ${root}/${s}/${sampleName}.findbest | head -1 | awk '{print $1}'` 
	echo "The best found: $best"
	#rm -r `ls | grep -v ${best}` 
	#find . -name "*gz" -exec mv {} ./ \;
	#find . -name "*fasta" -exec mv {} ./ \;
	#rm -r ${best}
	done
#else
#    for i in `ls $allref`; do cp $i ./; done
#    for i in *fasta; do
#        (cd ${root}
#        mkdir ${i%.fasta}
#        cp ${i} ${i%.fasta}
#        cp *fastq ${i%.fasta}
#        echo "working on $sampleName $i"
#        cd ${i%.fasta}; findbest) &
#    let count+=1
#    [[ $((count%10)) -eq 0 ]] && wait
#    done
#fi
wait


# created 2015-08-10 stuber
