#!/bin/sh

root=`pwd`

#PATHs
picardPath='/usr/local/bin/picard-tools-1.117/'
GATKPath='/usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar'
pythonGetFasta="/home/tstuber/workspace/stuber/python_scripts/GetFASTAbyGI.py"

# NCBI downloaded reference location
mydb="/data/mydb"

# flag -m will email just "M"e
# flag -b will turn off muliple for starts for "B"ug finding
# flag -k will run Kraken
bflag=
mflag=
kflag=
while getopts 'bmk' OPTION; do
    case $OPTION in
        b) bflag=1
        ;;
        m) mflag=1
        ;;
        k) kflag=1
        ;;
        ?) echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))

# This must be below the getopts
argUsed=`echo $1 | tr '[:lower:]' '[:upper:]'`

#######################################################################################
#|||||||||||||||||||||||||||||| Environment Controls ||||||||||||||||||||||||||||||||||
#######################################################################################

if [[ $1 == sivall ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/std/"
    targetref=/home/shared/databases/viral/ai/sivall/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/ai/sivall/newfiles"
    echo "vcftofasta.sh ran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == aiall ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/std/"
    targetref=/home/shared/databases/viral/ai/aiall/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/ai/aiall/newfiles"
    echo "vcftofasta.sh ran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == h5n2 ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/std/"
    pingyrdb=yes #(yes or no) Do you want to BLAST pintail gyrfalcon database
    targetref=/home/shared/databases/viral/ai/h5n2/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/ai/h5n2/newfiles"
    echo "vcftofasta.sh rran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == h5n8 ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/std/"
    pingyrdb=yes #(yes or no) Do you want to BLAST pintail gyrfalcon database
    targetref=/home/shared/databases/viral/ai/h5n8/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/ai/h5n8/newfiles"
    echo "vcftofasta.sh rran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == secd ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/stdPlusSECD"
    targetref=/home/shared/databases/viral/secd/*fasta
    bioinfoVCF="/bioinfo11/TStuber/Results/viruses/secd/newfiles"
    echo "vcftofasta.sh ran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == reo ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/std/"
    targetref=/home/shared/databases/viral/reo/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/reo/newfiles"
    echo "vcftofasta.sh ran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == vsv ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/std/"
    targetref=/home/shared/databases/viral/vsv/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/vsv/newfiles"
    echo "vcftofasta.sh ran targeting $1"
    echo "Script vcftofasta.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

else
    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: aiall, sivall, h5n2, h5n8, secd, reo, vsv"
    echo ""
    echo "Set optional flags"
    echo -e '   flag -m will email just "M"e'
    echo -e '   flag -b will run for de"B"ugging'
    echo -e '   flag -k wll run "K"raken'
    echo ""
    echo "Example: [prompt]$ id_virus_reads.sh -mbk aiall"
    echo ""
    exit 1

fi

echo "Core number being used: $NR_CPUS"
sleep 4

mkdir originalreads
cp *fastq* originalreads

# Unzip files if needed, else put std error to null
find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

#Set sample name- sample name is the name of the fastq files minus any identifying information from the sequencer
sampleName=`ls *.fastq | head -1 | sed 's/_.*//' | sed 's/\..*//'`

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

echo "" > $root/${sampleName}.summaryfile
echo "" > $root/${sampleName}.detailfile
summaryfile="$root/${sampleName}.summaryfile"
detailfile="$root/${sampleName}.detailfile"
emailbody=${root}/emailfile
echo "Sample: ${sampleName}" >> $summaryfile
echo "Sample: ${sampleName}" >> $detailfile
echo "Sample: ${sampleName}" >> $emailbody
echo "Reference_Set: $argUsed" >> $summaryfile
echo "Reference_Set: $argUsed" >> $detailfile
echo "Reference_Set: $argUsed" >> $emailbody

#######################################################################################
#######################################################################################
# /////////////////////////////// SCRIPT STARTING POINT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#######################################################################################
#######################################################################################

echo "*** PARAMETER SETTINGS ***"
echo "Path to genotyping codes: $genotypingcodes"
echo "Path to references: $targetref"
echo "Path to upload files: $bioinfoVCF"
echo "email_list: $email_list"
echo "sample type: $sampleType"

echo "Core number being used: $NR_CPUS"
sleep 4

#Calculate files sizes
forFileSize=`ls -lh $forReads | awk '{print $5}'`
if [[ $sampleType == "paired" ]]; then
    revFileSize=`ls -lh $revReads | awk '{print $5}'`
fi

#Calculate count of reads in each file
forCount=`grep -c '^+$' $forReads`
echo "R1 read count: $forCount"
if [[ $sampleType == "paired" ]]; then
    revCount=`grep -c '^+$' $revReads`
fi

if [[ $sampleType == "paired" ]]; then
    echo "R1 file size: ${forFileSize}, read count: $forCount" >> $summaryfile
    echo "R2 file size: ${revFileSize}, read count: $revCount" >> $summaryfile

    echo "R1 file size: ${forFileSize}, read count: $forCount" >> ${emailbody}
    echo "R2 file size: ${revFileSize}, read count: $revCount" >> ${emailbody}
else
    echo "Single fastq file size: ${forFileSize}, read count: $forCount" >> $summaryfile
    echo "Single fastq file size: ${forFileSize}, read count: $forCount" >> $emailbody
fi

#######################################################################################
#|||||||||||||||||||||||||||||||||||| Kraken ||||||||||||||||||||||||||||||||||||||||||
#######################################################################################

if [ "$kflag" ]; then
    echo "Kraken database selected is: $krakenDatabase"
    date
    echo "*** Kraken is finding reads"
    #Run Kraken
    if [[ $sampleType == "paired" ]]; then
    kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired *fastq* > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-kraken_report.txt
    else
    kraken --db ${krakenDatabase} --threads ${NR_CPUS}  $forReads > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-kraken_report.txt
    fi
    date
    echo "*** Krona transforming Kraken output to graph"

    # Run Krona
    cut -f2,3 $sampleName-output.txt > $sampleName-kronaInput.txt;
    /usr/local/bin/ktImportTaxonomy $sampleName-kronaInput.txt;
    mv taxonomy.krona.html $sampleName-Krona_identification_graphic.html;
    mv taxonomy.krona.html.files $sampleName-taxonomy.krona.html.files

    # Set variables and paths
    output=`ls *-output.txt`
    report=`ls *kraken_report.txt`

    printf "%s, %s file size, %'.0f reads\n" ${forFile} ${forFileSize} ${forCount}

    if [[ $sampleType == "paired" ]]; then
        printf "%s, %s file size, %'.0f reads\n" ${revFile} ${revFileSize} ${revCount}
    fi

    declare -i x=${forCount}
    declare -i y=${revCount}

    echo "" | awk -v x=$x -v y=$y '{printf "Total single end read count: %'\''d\n", x+y}'

    #Section of results summary that calculates number of reads per type of organism (ex: ssRNA virus)
    echo "Summary of Kraken Findings"
    cRead=`grep -c "^C" $output`
    uRead=`grep -c "^U" $output`
    virusreport=`awk ' $5 == "10239" {print $2}' $report`
    let allReads=cRead+uRead
    echo "allReads: $allReads"

    if [ -z $virusreport ]; then
        virusreport="zero"
    fi
    declare -i v=${virusreport}
    declare -i z=${allReads}
    echo "v is $v"
    echo "z is $z"

    pvRead=`awk -v v=$v -v z=$z 'BEGIN { print (v / z)*100 }'`

    echo "`printf "%'.0f\n" ${virusreport}` virus reads --> ${pvRead}% of total reads" >> $summaryfile
    echo "`printf "%'.0f\n" ${virusreport}` virus reads --> ${pvRead}% of total reads" >> $emailbody
else
    echo "Kraken not ran"
fi

#######################################################################################
#||||||||||||||||||||||||||||| Function to align reads ||||||||||||||||||||||||||||||||
#######################################################################################
function alignreads () {

echo ""
echo " ********************** ALIGNMENT INTERATION 1 STARTED ********************** "
echo ""

ref=`ls | grep .fasta`
echo "Reference Input:  $ref"
orgref=${ref%.fasta}
refname=${ref%.fasta}
##
bwa index $ref
samtools faidx $ref
java -jar ${picardPath}/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -jar ${picardPath}/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict
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
rm ${refname}.sam
rm ${refname}.raw.bam
samtools view -h -b -F4 ${refname}.sorted.bam > ./$refname.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picardPath}/MarkDuplicates.jar INPUT=${refname}.mappedReads.bam OUTPUT=$refname.dup.bam METRICS_FILE=$refname.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $refname.dup.bam"
samtools index $refname.dup.bam

echo "*** Calling VCFs with UnifiedGenotyper"
java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I $refname.dup.bam -o ${refname}.UG.vcf -nct 2

if [ -s ${refname}.UG.vcf ]; then
    echo "${g}.UG.vcf present continueing script"
else
    sleep 60
    echo "${refname}.UG.vcf is missing, try making again"
    java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${dupbam} -o ${refname}.UG.vcf -nct 2
    if [ -s ${refname}.UG.vcf ]; then
        echo "${refname}.UG.vcf present continueing script"
    else
        echo "" >> ${emailbody}
        echo "${refname}.UG.vcf failed to make, Further analysis required" >> ${emailbody}
        email_list="Tod.P.Stuber@aphis.usda.gov"
        echo "${refname}.UG.vcf failed to make, Further analysis required" | mutt -s "Sample $sampleName Reference_Set $argUsed" -- $email_list
        exit 1
    fi
fi

# make reference guided contig
java -jar ${GATKPath} -T FastaAlternateReferenceMaker -R $ref -o ${refname}.readreference.fasta -V ${refname}.UG.vcf

echo ">${refname}" > ${refname}.readreference.temp; awk ' $8 ~ /^AN=2/ {print $4} ' ${refname}.UG.vcf | tr -d [:space:] >> ${refname}.readreference.temp; echo "" >> ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.readreference.fasta

echo "short BLAST"
blastn -query ${refname}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${refname}-readreference-max1-nt-id.txt -max_target_seqs 1 -outfmt "6 saccver"

if [ -s ${refname}-readreference-max1-nt-id.txt ]; then
    echo "Something was found"
else
    echo "No matches"
    exit 1
fi

rm *dict
rm *fasta*
rm *coveragefile
rm *vcf*
rm *bam*

#######################################################################################
# End of initial iteration

echo "$refname"
echo " ********************** ALIGNMENT INTERATION 2 STARTED ********************** "
echo ""
alignagain

echo "$refname"
echo " ********************** ALIGNMENT INTERATION 3 STARTED ********************** "
echo ""
alignagain

# Start of last iteration
#######################################################################################

echo "$refname"
echo " ********************** ALIGNMENT INTERATION 4 STARTED ********************** "
echo ""

acc=`head -1 ${refname}-readreference-max1-nt-id.txt | sed 's/\..*//'`
writeFile="${acc}.fasta"
echo "acc $acc"
echo "writeFile $writeFile"
pwd

python -u ${pythonGetFasta} $acc $writeFile
wait

ls ${mydb} >  ${mydb}/list.txt
p=`grep "${acc}" ${mydb}/list.txt`

if [[ -z "$p" ]]; then
    echo "Downloading from NCBI"
    echo "This is the file to write to:  $writeFile"
    echo "Running python script to grab fasta file"
    grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
    if [ -s $writeFile ]; then
        echo "Downloaded from NCBI, Good to continue."
    else
        echo "Try downloading again"
        sleep 20
        echo "Running python script to grab fasta file"
        python -u ${pythonGetFasta} $acc $writeFile
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        sleep 5
        if [ ! -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 120
            echo "Running python script to grab fasta file"
            python -u ${pythonGetFasta} $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 320
                echo "Running python script to grab fasta file"
                python -u ${pythonGetFasta} $acc $writeFile
                grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                sleep 5
                if [ ! -s $writeFile ]; then
                    read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
                fi
            fi
        fi
    fi
    cp $writeFile ${mydb}
else
    echo "File is local"
    cp ${mydb}/${p} ./
    writeFile=${p}
fi

# Grab reads and reference and place them in variables
ref=${writeFile}
writefilelist=${root}/writelist
grep ">" $writeFile >> $writefilelist
export $writefilelist

echo "Reference Input:  $ref"

refname=${ref%.fasta}
echo "*** refname $refname"

lastalign

###
# Place zero coverage positins to hapreadyAll.vcf

echo "${orgref}-${refname}.hapreadyAll.vcf"
echo "${orgref}-${refname}.UG.vcf"

grep '^#' ${orgref}-${refname}.hapreadyAll.vcf > header
grep -v '^#' ${orgref}-${refname}.hapreadyAll.vcf > snps

awk 'BEGIN{OFS="\t"} $10 == "./." {print $0}' ${orgref}-${refname}.UG.vcf > zerocoverage

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "N", $6, $7, $8, "GT", "1"}' zerocoverage > zeroformated

cat snps zeroformated | sort -nk2,2 > body

cat header body > ${orgref}-${refname}.hapreadyAll.vcf
###

# make reference guided contig
java -jar ${GATKPath} -T FastaAlternateReferenceMaker -R $ref -o ${orgref}-${refname}.reference_guided.fasta -V ${orgref}-${refname}.hapreadyAll.vcf -IUPAC ${orgref}-${refname}

if [ $sampleType == "paired" ]; then
    java -jar ${picardPath}/SamToFastq.jar INPUT=${orgref}-${refname}.dup.bam FASTQ=./${orgref}-${refname}-mapped_R1.fastq SECOND_END_FASTQ=./${orgref}-${refname}-mapped_R2.fastq
else
    java -jar ${picardPath}/SamToFastq.jar INPUT=${orgref}-${refname}.dup.bam FASTQ=./${orgref}-${refname}-mapped_R1.fastq
fi

if [ $sampleType == "paired" ]; then
    mapCount=`cat ${orgref}-${refname}-mapped_R1.fastq ${orgref}-${refname}-mapped_R2.fastq | grep -c "^+$"`
else
    mapCount=`cat ${orgref}-${refname}-mapped_R1.fastq | grep -c "^+$"`
fi
echo "mapCount: $mapCount"

#Length of reference
countNTs=`grep -v ">" $ref | wc | awk '{print $3}'`

#Number of nucleotides in reference with coverage
echo "*** Bamtools is getting coverage..."
bamtools coverage -in ${orgref}-${refname}.dup.bam | awk -v x=${orgref}-${refname} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${orgref}-${refname}-coveragefile

covCount=`awk '{ if ($3 != 0) count++ } END { print count }' ${orgref}-${refname}-coveragefile`
echo "covCount $covCount"

declare -i x=${covCount}
declare -i y=${countNTs}

#Percent of reference with coverage
perc=`awk -v x=$x -v y=$y 'BEGIN { print(x/y)*100}'`
echo "perc: $perc"

#Average depth of coverage for the reference
depthcov=`awk '{sum+=$3} END { print sum/NR}' ${orgref}-${refname}-coveragefile`
echo "depthcov $depthcov"

LC_NUMERIC=en_US
printf "%-30s %'10d %11.2f%% %'10dX\n" ${orgref}-${refname} $mapCount $perc $depthcov >> ${summaryfile}.pre

echo "`printf "%'.0f\n" ${mapCount}` reads aligned to ${orgref}-${refname}"
echo "${perc}% genome coverage, $depthcov"

awk -v x=${orgref}-${refname} 'BEGIN {OFS="\t"} {print x, $2, $3}' ${orgref}-${refname}-coveragefile | grep -v "segment_all" > ${orgref}-${refname}-samplecoveragefile

rm *sorted.bam*
rm *dup.bam*
rm *mappedReads.bam
rm *fasta.amb
rm *fasta.ann
rm *fasta.bwt
rm *fasta.pac
rm *fasta.sa
rm *fastq*
#rm body
#rm header
#rm ${g}.cutlist
mkdir alignment
mv *fastq alignment
mv *dict alignment
mv *fasta* alignment
mv *bam* alignment
mv *vcf* alignment
mv *pdf alignment
mv *coveragefile alignment

}


#######################################################################################
function lastalign () {
# Align Last Time
bwa index $ref
samtools faidx $ref
java -jar ${picardPath}/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -jar ${picardPath}CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${n}-${refname}.dict
    if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

#${orgref}-${refname}.sam

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allow
if [[ $sampleType == "paired" ]]; then
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads $revReads > ${orgref}-${refname}.sam
else
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads > ${orgref}-${refname}.sam
fi

samtools view -bh -F4 -T $ref ${orgref}-${refname}.sam > ${orgref}-${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${orgref}-${refname}.raw.bam ${orgref}-${refname}.sorted
echo "****Indexing Bam"
samtools index ${orgref}-${refname}.sorted.bam
rm ${orgref}-${refname}.sam
rm ${orgref}-${refname}.raw.bam
samtools view -h -b -F4 ${orgref}-${refname}.sorted.bam > ./${orgref}-${refname}.mappedReads.bam
samtools index ./${orgref}-${refname}.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picardPath}/MarkDuplicates.jar INPUT=${orgref}-${refname}.mappedReads.bam OUTPUT=${orgref}-${refname}.dup.bam METRICS_FILE=${orgref}-${refname}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${orgref}-${refname}.dup.bam"
samtools index ${orgref}-${refname}.dup.bam
java -jar ${GATKPath} -T ClipReads -R $ref -I ${orgref}-${refname}.dup.bam -o ${orgref}-${refname}.downsample.bam -filterNoBases -dcov 10
samtools index ${orgref}-${refname}.downsample.bam

rm *UG.vcf

#bam prepared now onto variant calling

java -Xmx4g -jar ${GATKPath} -R $ref -T HaplotypeCaller -ploidy 1 -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.hapreadyAll.vcf -bamout ${orgref}-${refname}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
java -Xmx4g -jar ${igvtools} index ${orgref}-${refname}.hapreadyAll.vcf

if [ -s ${orgref}-${refname}.hapreadyAll.vcf ]; then
    echo "${orgref}-${refname}.hapreadyAll.vcf present continueing script"
else
    sleep 60
    echo "${orgref}-${refname}.hapreadyAll.vcf is missing, try making again"
    java -Xmx4g -jar ${GATKPath} -R $ref -T HaplotypeCaller -ploidy 1 -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.hapreadyAll.vcf -bamout ${orgref}-${refname}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
    java -Xmx4g -jar ${igvtools} index ${orgref}-${refname}.hapreadyAll.vcf
    if [ -s ${orgref}-${refname}.hapreadyAll.vcf ]; then
        echo "${orgref}-${refname}.hapreadyAll.vcf present continueing script"
    else
        echo "" >> ${emailbody}
        echo "${orgref}-${refname}.hapreadyAll.vcf failed to make, Further analysis required" >> ${emailbody}
        email_list="Tod.P.Stuber@aphis.usda.gov"
        echo "${orgref}-${refname}.hapreadyAll.vcf failed to make, Further analysis required" | mutt -s "Sample $sampleName Reference_Set $argUsed" -- $email_list
        exit 1
    fi
fi

echo "*** Calling VCFs with UnifiedGenotyper"
java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.UG.vcf -nct 2

if [ -s ${orgref}-${refname}.UG.vcf ]; then
    echo "${orgref}-${refname}.UG.vcf present continueing script"
else
    sleep 60
    echo "${orgref}-${refname}.UG.vcf is missing, try making again"
    java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.UG.vcf -nct 2
    if [ -s ${orgref}-${refname}.UG.vcf ]; then
        echo "${orgref}-${refname}.UG.vcf present continueing script"
    else
        echo "" >> ${emailbody}
        echo "${orgref}-${refname}.UG.vcf failed to make, Further analysis required" >> ${emailbody}
        email_list="Tod.P.Stuber@aphis.usda.gov"
        echo "${orgref}-${refname}.UG.vcf failed to make, Further analysis required" | mutt -s "Sample $sampleName Reference_Set $argUsed" -- $email_list
        exit 1
    fi
fi

}

#######################################################################################
#|||||||||||||||||||||||||||||| Interative Alignment ||||||||||||||||||||||||||||||||||
#######################################################################################

function alignagain () {

acc=`head -1 ${refname}-readreference-max1-nt-id.txt | sed 's/\..*//'`
writeFile="${acc}.fasta"
echo "acc $acc"
echo "writeFile $writeFile"
pwd

python -u ${pythonGetFasta} $acc $writeFile
wait

ls ${mydb} >  ${mydb}/list.txt
p=`grep "${acc}" ${mydb}/list.txt`

if [[ -z "$p" ]]; then
    echo "Downloading from NCBI"
    echo "This is the file to write to:  $writeFile"
    echo "Running python script to grab fasta file"
    grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
    if [ -s $writeFile ]; then
        echo "Downloaded from NCBI, Good to continue."
    else
        echo "Try downloading again"
        sleep 20
        echo "Running python script to grab fasta file"
        python -u ${pythonGetFasta} $acc $writeFile
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        sleep 5
        if [ ! -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 120
            echo "Running python script to grab fasta file"
            python -u ${pythonGetFasta} $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 320
                echo "Running python script to grab fasta file"
                python -u ${pythonGetFasta} $acc $writeFile
                grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                sleep 5
                if [ ! -s $writeFile ]; then
                    read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
                fi
            fi
        fi
    fi
    cp $writeFile ${mydb}
else
    echo "File is local"
    cp ${mydb}/${p} ./
    writeFile=${p}
fi

#grep ">" $writeFile >> $summaryfile
refheader=`grep ">" $writeFile`

# Grab reads and reference and place them in variables
ref=${writeFile}

echo "Reference Input:  $ref"
refname=${ref%.fasta}

bwa index $ref
samtools faidx $ref
java -jar ${picardPath}/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -jar ${picardPath}CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict
    if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allowed
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
rm ${refname}.sam
rm ${refname}.raw.bam

echo "*** Calling VCFs with UnifiedGenotyper"
java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${refname}.sorted.bam -o ${refname}.UG.vcf -nct 2

if [ -s ${refname}.UG.vcf ]; then
    echo "${refname}.UG.vcf present continueing script"
else
    sleep 60
    echo "${refname}.UG.vcf is missing, try making again"
    java -jar ${GATKPath} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${refname}.sorted.bam -o ${refname}.UG.vcf -nct 2
    if [ -s ${refname}.UG.vcf ]; then
        echo "${refname}.UG.vcf present continueing script"
    else
        echo "" >> ${emailbody}
        echo "${refname}.UG.vcf failed to make, Further analysis required" >> ${emailbody}
        email_list="Tod.P.Stuber@aphis.usda.gov"
        echo "${refname}.UG.vcf failed to make, Further analysis required" | mutt -s "Sample $sampleName Reference_Set $argUsed" -- $email_list
        exit 1
    fi
fi

# make reference guided contig
java -jar ${GATKPath} -T FastaAlternateReferenceMaker -R $ref -o ${refname}.readreference.fasta -V ${refname}.UG.vcf

sed 's/NN//g' < ${refname}.readreference.fasta > ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.readreference.fasta

echo "short BLAST"
blastn -query ${refname}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${refname}-readreference-max1-nt-id.txt -max_target_seqs 1 -outfmt "6 saccver"

rm *dict
rm *fasta*
rm *coveragefile
rm *vcf*
rm *bam*

}

#######################################################################################
#|||||||||||||||||||||| Function to make R graph of coverage ||||||||||||||||||||||||||
#######################################################################################
function plotR () {
cat > ./plotR.r << EOL
#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)

arg <- commandArgs(trailingOnly=TRUE)

data <- read.csv(arg[1], header=FALSE, sep="\t")
names(data) <- c("species", "position", "coverage")

pdf("myplot.pdf", width=20, height=4)

ggplot(data, aes(x=position, y=log(coverage), colour=species, group=species)) + geom_point(size=2.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()
EOL

chmod 755 ./plotR.r

./plotR.r $1 $2
rm ./plotR.r
}

#######################################################################################
#||||||||||||||||||||||||||||||||| Reference Setup ||||||||||||||||||||||||||||||||||||
#######################################################################################

echo "Set up references"
# Reference is set in Environment Controls
# target will only call "fasta" files

if [ "$bflag" ]; then
    echo ""
    echo " *** B FLAG ON, BUG FINDING MODE, SINGLE SAMPE PROCESSING *** "
    echo ""
    for i in `ls $targetref`; do cp $i ./; done
    for i in *fasta; do
        cd ${root}
        mkdir ${i%.fasta}
        cp ${i} ${i%.fasta}
        cp *fastq ${i%.fasta}
        echo "working on $sampleName $i"
        cd ${i%.fasta}; alignreads
    done
else
    for i in `ls $targetref`; do cp $i ./; done
    for i in *fasta; do
        (cd ${root}
        mkdir ${i%.fasta}
        cp ${i} ${i%.fasta}
        cp *fastq ${i%.fasta}
        echo "working on $sampleName $i"
        cd ${i%.fasta}; alignreads) &
    let count+=1
    [[ $((count%10)) -eq 0 ]] && wait
    done
fi
wait

#######################################################################################
#|||||||||||||||||||||||||||||||||||| Finish Up |||||||||||||||||||||||||||||||||||||||
#######################################################################################

wait
sleep 5
wait
sleep 5
echo "" >> ${summaryfile}
echo "Alignment stats (reference guided):" >> ${summaryfile}
printf "%-30s %10s %11s %10s\n" "reference used" "read count" "percent cov" "ave depth" >> ${summaryfile}
sort -k1,1 < ${summaryfile}.pre >> ${summaryfile}

echo "" >> ${emailbody}
echo "Alignment stats (reference guided):" >> ${emailbody}
printf "%-30s %10s %11s %10s\n" "reference used" "read count" "percent cov" "ave depth" >> ${emailbody}
sort -k1,1 < ${summaryfile}.pre >> ${emailbody}

rm ${summaryfile}.pre
#####################

for i in `find . -name "*samplecoveragefile"`; do
    cat $i >> allsamplecoveragefile
done

plotR allsamplecoveragefile $sampleName
mv myplot.pdf ${sampleName}.assembly_graph.pdf

cd $root
pwd

mkdir ${sampleName}-reference_guided_assemblies
for i in `find . -name "*reference_guided.fasta"`; do
    cp $i ${sampleName}-reference_guided_assemblies
done

cd ${root}/${sampleName}-reference_guided_assemblies
pwd

for i in `ls *fasta | sort -k1,1`; do
    echo ">${i%.reference_guided.fasta}" >> ${sampleName}.consensus.reads.fasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | awk 'length > x { x = length; y = $0 } END { print y }' >> ${sampleName}.consensus.reads.fasta
done

#######################################################################################
#|||||||||||||||||||||||||||||| Reference Set Alignment |||||||||||||||||||||||||||||||
#######################################################################################

assemblyfolder=`pwd`

mkdir ${root}/igv_alignment
cp ${sampleName}.consensus.reads.fasta ${root}/igv_alignment
cd ${root}/igv_alignment

echo "forward read: $forReads"
echo "reverse read: $revReads"

ref=`ls | grep .fasta`
echo "Reference Input:  $ref"
refname=${ref%.fasta}
orgref="igvalignment"

#############################

bwa index $ref
samtools faidx $ref
java -jar ${picardPath}/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
echo "Index and dict are present, continue script"
else
sleep 5
echo "Either index or dict for reference is missing, try making again"
samtools faidx $ref
java -jar ${picardPath}CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${n}-${refname}.dict
if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
fi
fi

#${orgref}-${refname}.sam

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allow
if [[ $sampleType == "paired" ]]; then
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads $revReads > ${orgref}-${refname}.sam
else
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads > ${orgref}-${refname}.sam
fi

samtools view -bh -F4 -T $ref ${orgref}-${refname}.sam > ${orgref}-${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${orgref}-${refname}.raw.bam ${orgref}-${refname}.sorted
echo "****Indexing Bam"
samtools index ${orgref}-${refname}.sorted.bam
rm ${orgref}-${refname}.sam
rm ${orgref}-${refname}.raw.bam
samtools view -h -b -F4 ${orgref}-${refname}.sorted.bam > ./${orgref}-${refname}.mappedReads.bam
samtools index ./${orgref}-${refname}.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picardPath}/MarkDuplicates.jar INPUT=${orgref}-${refname}.mappedReads.bam OUTPUT=${orgref}-${refname}.dup.bam METRICS_FILE=${orgref}-${refname}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${orgref}-${refname}.dup.bam"
samtools index ${orgref}-${refname}.dup.bam
java -jar ${GATKPath} -T ClipReads -R $ref -I ${orgref}-${refname}.dup.bam -o ${orgref}-${refname}.downsample.bam -filterNoBases -dcov 10
samtools index ${orgref}-${refname}.downsample.bam

#############################

rm *sorted.bam*
rm *dup.bam*
rm *mappedReads.bam
rm *fasta.amb
rm *fasta.ann
rm *fasta.bwt
rm *fasta.pac
rm *fasta.sa
rm *fastq*
rm *dict
rm *mappedReads*

# Go back to the assemble folder
cd $assemblyfolder

###########

contigcount=`grep -c ">" ${sampleName}.consensus.reads.fasta`
sed 's/NNN//g' ${sampleName}.consensus.reads.fasta > ${sampleName}.consensusnoN.reads.fasta

echo "nt BLAST $contigcount contigs..."
blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${sampleName}-consensus-max1-nt.txt -max_target_seqs 1 -outfmt "6 qseqid qlen slen pident mismatch evalue bitscore stitle saccver"
echo "" >> ${summaryfile}
echo "--------------------------------------------------" >> ${summaryfile}
echo "*** NT Database ***" >> ${summaryfile}
awk 'BEGIN{printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s\n", "query ID", "qlength", "slength", "% id", "mis", "evalue", "bscore", "Description"} {printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' ${sampleName}-consensus-max1-nt.txt >> ${summaryfile}
echo "" >> ${summaryfile}

if [ $pingyrdb == yes ]; then
    echo "pintail-gyrfalcon BLAST $contigcount contigs..."
    blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/blastdb-pintail-gyrfalcon/pintail-gyrfalcon.fsa -num_threads 20 -out ${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt -outfmt 1
    blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/blastdb-pintail-gyrfalcon/pintail-gyrfalcon.fsa -num_threads 20 -out ${sampleName}-consensus-fmt6-pintail-gyrfalcon.txt -outfmt "6 qseqid qlen slen pident mismatch evalue bitscore stitle saccver"
    echo "--------------------------------------------------" >> ${summaryfile}
    echo "*** Pintail/Gyrfalcon Database ***" >> ${summaryfile}
    awk 'BEGIN{printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s\n", "query ID", "qlength", "slength", "% id", "mis", "evalue", "bscore", "Description"} {printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' ${sampleName}-consensus-fmt6-pintail-gyrfalcon.txt >> ${summaryfile}
fi

###########################

mkdir ${root}/share_folder
cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/share_folder

############################
cd ${root}/share_folder

#Format fasta headers for NCBI

if [[ ${genotypingcodes} == "NEED TO SET" ]]; then
    echo "genotyping codes not given"
    cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/share_folder/${sampleName}-submissionfile.fasta
else
    echo "sampleName: $sampleName"

    grep "$sampleName" $genotypingcodes | head -1 > ${sampleName}.information
    cat ${sampleName}.information

    if [ -s ${sampleName}.information ]; then
        echo "file exists and is greater than zero, continue"

        #column 1: sample
        sample=`awk 'BEGIN{FS="\t"}{print $1}' ${sampleName}.information`
        echo "sample $sample"

        #column 2: species
        species=`awk 'BEGIN{FS="\t"}{print $2}' ${sampleName}.information`
        echo "species $species"

        #column 3: state
        state=`awk 'BEGIN{FS="\t"}{print $3}' ${sampleName}.information`
        echo "state $state"

        syear=`echo "$sample" | sed 's/-.*//'`
        echo "syear $syear"
        sampleyear=`echo "20${syear}"`
        echo "sampleyear $sampleyear"

        sed 's/>.*seg.*1_.*/>Seq1/' ${sampleName}.consensus.reads.fasta | sed 's/>.*seg.*2_.*/>Seq2/' | sed 's/>.*seg.*3_.*/>Seq3/' | sed 's/>.*seg.*4_.*/>Seq4/' | sed 's/>.*seg.*5_.*/>Seq5/' | sed 's/>.*seg.*6_.*/>Seq6/' | sed 's/>.*seg.*7_.*/>Seq7/' | sed 's/>.*seg.*8_.*/>Seq8/' > ${sampleName}.temp

        noyear=`echo $sample | sed -e "s/$syear-//"`
        echo "noyear $noyear"

        # Create "here-document"
cat >./param.txt <<EOL

>Seq1 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 1, polymerase PB2 (PB2) gene, complete cds.
>Seq2 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 2, polymerase PB1 (PB1) gene, complete cds.
>Seq3 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 3, polymerase PA (PA) gene, complete cds.
>Seq4 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 4, hemagglutinin (HA) gene, complete cds.
>Seq5 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 5, nucleoprotein (NP) gene, complete cds.
>Seq6 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 6, neuraminidase (NA) gene, complete cds.
>Seq7 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 7, matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds.
>Seq8 [organism=Influenza A virus](A/${species}/${state}/${noyear}/${sampleyear}(${argUsed})) segment 8, non-structural protein NS1 and non-structural protein NS2 (NS) gene, complete cds.

EOL

        awk 'NR==FNR{a[$1]=$0;next} ($1 in a){ print a[$1]; next}1' param.txt ${sampleName}.temp > ${sampleName}-submissionfile.fasta

    else
        echo "metadata not available"
        cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/share_folder/${sampleName}-submissionfile.fasta
    fi

rm *temp
rm *information
rm param.txt

fi

###########################
cd $root
#rm *fastq*
echo "" >> ${summaryfile}
echo "" >> ${summaryfile}
echo "Files copied to: $bioinfoVCF" >> ${summaryfile}
echo "" >> ${summaryfile}

rm *headers
rm *kronaInput.txt
#rm allsamplecoveragefile
rm -r *taxonomy.krona.html.files
rm ${sampleName}.detailfile

mkdir kraken
mv *output.txt kraken
mv *report.txt kraken
cp $0 ${root}
echo "******* $LINENO, $PWD"
fileName=`basename $0`
cp ${root}/share_folder/${sampleName}-submissionfile.fasta ${root}
cp ${sampleName}-reference_guided_assemblies/${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt ${root}

enscript ${summaryfile} -B -j -r -f "Courier5" -o - | ps2pdf - ${sampleName}-report.pdf

if [ -e ${sampleName}-report.pdf ]; then
	ls ${sampleName}-report.pdf > emailfiles
fi

if [ -e ${sampleName}.assembly_graph.pdf ]; then
	ls ${sampleName}.assembly_graph.pdf >> emailfiles
fi

if [ -e ${root}/${sampleName}-submissionfile.fasta ]; then
	ls ${root}/${sampleName}-submissionfile.fasta >> emailfiles
fi

if [ -e ${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt ]; then
	ls ${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt >> emailfiles
fi

if [ -e $sampleName-Krona_identification_graphic.html ]; then 
	ls $sampleName-Krona_identification_graphic.html >> emailfiles
fi

if [[ $sampleType == "paired" ]]; then
	echo "paried data, not checking for C insert"
else
	if [ $pingyrdb == yes ]; then
		noc=`egrep -c "GAGTTGACATAAACCAGGCCACGC|GCGTGGCCTGGTTTATGTCAACTC" $forReads`
		cinsert=`egrep -c "GAGTTGACATAAACCCAGGCCACGC|GCGTGGCCTGGGTTTATGTCAACTC" $forReads`
        insert1=`egrep -c "GAGTTGACATAAA[AGT]CCAGGCCACGC|GCGTGGCCTGG[ACT]TTTATGTCAACTC" $forReads`
        insert2=`egrep -c "GAGTTGACATAAAC[AGT]CAGGCCACGC|GCGTGGCCTG[ACT]GTTTATGTCAACTC" $forReads`
        insert3=`egrep -c "GAGTTGACATAAACC[AGT]AGGCCACGC|GCGTGGCCT[ACT]GGTTTATGTCAACTC" $forReads`
		echo "" >> ${emailbody}
		echo "(-CC) No insert count: $noc" >> ${emailbody}
		echo "(CCC) C insert count: $cinsert" >> ${emailbody}
        echo "([AGT]CC) Non-C insert count at position 1: $insert1" >> ${emailbody}
        echo "(C[AGT]C) Non-C insert count at position 2: $insert2" >> ${emailbody}
        echo "(CC[AGT]) Non-C insert count at position 3: $insert3" >> ${emailbody}
		echo "" >> ${emailbody}
	else
		echo "pingyrdb not being referenced, therefore not checking for C insert"
	fi
fi
rm *fastq*

Cleanup
rm -r `ls | egrep -v "emailfile|emailfiles|$0|igv_alignment|originalreads|summaryfile|report.pdf|Krona_identification_graphic.html|-consensus-blast_alignment-pintail-gyrfalcon.txt|-submissionfile.fasta|assembly_graph.pdf"`

if [ "$mflag" ]; then
    email_list="tod.p.stuber@usda.gov"
    cat ${emailbody} | mutt -s "Sample: ${sampleName}, Reference_Set: $argUsed" -a `cat emailfiles` -- $email_list
    rm emailfiles
else
    cat ${emailbody} | mutt -s "Sample: ${sampleName}, Reference_Set: $argUsed" -a `cat emailfiles` -- $email_list

    rm ${emailbody}
    rm emailfiles
    echo "Copying to ${bioinfoVCF}"
    cp -r $PWD ${bioinfoVCF}

fi

echo "****************************** END ******************************"
pwd


# Created: 2015-06-12, tstuber
