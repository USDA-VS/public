#!/bin/sh

alias pause='read -p "$LINENO Enter"'

root=`pwd`
flu=no

#PATHs
picard='/usr/local/bin/picard-tools-1.141/picard.jar'
gatk='/usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar'
pythonGetFasta=`which GetFASTAbyGI.py`

# NCBI downloaded reference location
mydb="/data/mydb"
#Delete files in local database that may be empty
dzdo chmod 755 * ${mydb}/*
for i in ${mydb}/*; do
	if [ -s $i ]; then
        	echo "file $i exists"
	else
        	echo "file $i is empty and has been deleted"
		rm -f $i
	fi
done

idscriptrunsummary="/home/shared/idvirus_run_summary.txt"

# flag -m will email just "M"e
# flag -b will turn off muliple for starts for "B"ug finding
# flag -k will run Kraken
# flag -e flag used when running script from idemail.sh
bflag=
mflag=
kflag=
eflag=
while getopts 'bmke' OPTION; do
    case $OPTION in
        b) bflag=1
        ;;
        m) mflag=1
        ;;
        k) kflag=1
        ;;
        e) eflag=1
        ;;
        ?) echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))

# This must be below the getopts
argUsed=`echo $1 | tr '[:lower:]' '[:upper:]'`

pingyrdb=""

# By default

#####################################################

function getfluname () {

echo "using xlrd to get flu genotyping codes from MiSeq_Samples.xlsx"
date

cat >./excelcolumnextract.py <<EOL
#!/usr/bin/env python

import xlrd
from sys import argv

script, input = argv

wb = xlrd.open_workbook(input)

sh = wb.sheet_by_index(0)
for rownum in range(sh.nrows):
    row = sh.row_values(rownum)
    # use the join/map/str to output comma delimited list only
    print (', '.join(map(str, row)))

EOL

chmod 755 ./excelcolumnextract.py

./excelcolumnextract.py /bioinfo11/MKillian/MiSeq\ samples/MiSeq_Samples.xlsx | sed 's/, /,/g' | awk 'BEGIN{FS=","; OFS="\t"} {print $2, $4, $5}' | sed -e 's/[.*:()/\?]/_/g' -e 's/ /_/g' -e 's/_-/_/' -e 's/-_/_/' -e 's/__/_/g' -e 's/[_-]$//' > /scratch/report/flu_genotyping_codes.txt 

rm ./excelcolumnextract.py

}

#####################################################


#######################################################################################
#|||||||||||||||||||||||||||||| EnvironmentControls ||||||||||||||||||||||||||||||||||
#######################################################################################

if [[ $1 == flu ]]; then
    getfluname
    flu=yes
    genotypingcodes="/scratch/report/flu_genotyping_codes.txt"
    krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts"
    pingyrdb=yes #(yes or no) Do you want to BLAST pintail gyrfalcon database
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/ai/flu/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/influenza/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov" 

elif [[ $1 == secd ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/secd/*fasta
    bioinfoVCF="/bioinfo11/TStuber/Results/viruses/secd/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == reo ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/reo/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/reo/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == vsv ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/vsv/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/vsv/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == isav ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/isav/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/isav/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == bvd ]]; then
    genotypingcodes="NEED TO SET"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/bvd/*fasta
    bioinfoVCF="/bioinfo11/KBrien/newfiles"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Kaitlin.E.Brien@aphis.usda.gov"

elif [[ $1 == newcastle ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/newcastle/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/newcastle"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == prrsv ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/prrsv/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/prrsv"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == pox ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/pox/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/pox"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == herpes ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/herpes/*fasta
    bioinfoVCF="/bioinfo11/MKillian/Analysis/results/herpes"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov" #mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == rhabdo ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/rhabdoviruses/*fasta
    #bioinfoVCF="/bioinfo11/MKillian/Analysis/results/rhabdoviruses"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov" # Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == pentro ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/porcine_enterovirus/*fasta
    #bioinfoVCF="/bioinfo11/MKillian/Analysis/results/porcine_enterovirus"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov" # Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == tescho ]]; then
    genotypingcodes="/bioinfo11/MKillian/Analysis/results/genotypingcodes.txt"
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    targetref=/bioinfo11/MKillian/Analysis/script_dependents/teschovirus/*fasta
    #bioinfoVCF="/bioinfo11/MKillian/Analysis/results/porcine_enterovirus"
    echo "idvirus.sh ran targeting $1"
    echo "Script idvirus.sh ran targeting $1"
    email_list="tod.p.stuber@usda.gov" # Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

else
    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: flu, secd, reo, vsv, isav, bvd, newcastle, prrsv, herpes, rhabdo, tescho, pentro (porcine_entrovirus)"
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
find . -maxdepth 1 -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

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

    echo "*** Trimming"
    #######################
    #                     #
    #      Trimming       #
    #                     #
    #######################
    
    forReads=`ls | grep _R1`
    echo "Forward Reads to be trimmed: $forReads"
    
    revReads=`ls | grep _R2`
    echo "Reverse Reads to be trimmed:: $revReads"
    
    #Trim the reads with bbmap tool kit (bbduk plugin)
    #about twice as fast as trimmomatic
    
    strain=$(echo $revReads | sed 's/_.*//' | sed 's/\..*//')
    echo -e "Quality trimming sample "$strain""
    
        bbduk.sh -Xmx80g \
        in1="$forReads" \
        in2="$revReads" \
        ref="/usr/local/bin/bbmap/resources/nextera_flu.fa.gz" \
        ktrim=r k=23 mink=11 hdist=1 \
        qtrim=lr trimq=5 \
        minlen=36 \
        out1=trimmed_reads/${strain}_Trimmed_R1.fastq.gz \
        out2=trimmed_reads/${strain}_Trimmed_R2.fastq.gz \
        stats=trim_stats.txt \
        qchist=qc_by_base.txt \
        threads=auto \
        showspeed=f
    
    mv -v trimmed_reads/${strain}_Trimmed_R1.fastq.gz ./
    mv -v trimmed_reads/${strain}_Trimmed_R2.fastq.gz ./
    rm -r trimmed_reads
    rm "$forReads"
    rm "$revReads"
    
    forReads=`ls | grep _R1`
    echo "Forward Reads to be used after trimmed: $forReads"
    
    revReads=`ls | grep _R2`
    echo "Reverse Reads to be used after trimmed:: $revReads"

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
mytex="$root/${sampleName}.tex"
detailfile="$root/${sampleName}.detailfile"
emailbody=${root}/emailfile
echo "Sample: ${sampleName}" >> $summaryfile
echo "Sample: ${sampleName}" >> $detailfile
echo "Sample: ${sampleName}" >> $emailbody
echo "Reference_Set: $argUsed" >> $summaryfile
echo "Reference_Set: $argUsed" >> $detailfile
echo "Reference_Set: $argUsed" >> $emailbody

#here-document
#latex preamble
cat << EOL > $mytex
\documentclass[a4paper,10pt]{article}
\usepackage[margin=0.5in]{geometry}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{floatrow}
\usepackage{float}
\floatsetup[table]{capposition=top}
\usepackage{caption}
\captionsetup{labelformat=empty,justification=justified,singlelinecheck=false}
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}
\usepackage{lastpage}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}

\cfoot{Appendix --  page \thepage\ of \pageref{LastPage}}

\begin{document}

\includegraphics[scale=0.2]{/home/shared/usdalogo.png}


\today

\vspace{5mm}
\textbf{ } \\\\
\textbf{Identification Report:  ${sampleName}} \\\\ 
\textbf{XXXXXHNTYPEXXXXXXX} \\\\
\textbf{XXXXXSTRAINNAMEXXXXXXX} \\\\
EOL

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
forCount=`grep -c '^+$' $forReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "R1 read count: $forCount"
if [[ $sampleType == "paired" ]]; then
    revCount=`grep -c '^+$' $revReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
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

forsize=`ls -lh $forReads | awk '{print $5}'`
revsize=`ls -lh $revReads | awk '{print $5}'`
n=`echo $forFile | sed 's/[._].*//'`

echo "" >> ${mytex}.filestats
echo "\vspace{5mm}" >> ${mytex}.filestats
echo "" >> ${mytex}.filestats

echo "\begin{table}[H]" >> ${mytex}.filestats
echo "\begin{tabular}{ l | p{7cm} | p{7cm} }" >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "file name & $forFile & $revFile \\\\ " | sed "s/$n[._]//g" | sed 's/_/\\_/g' >> ${mytex}.filestats
echo "\hline \hline" >> ${mytex}.filestats
echo "read count & $forCount & $revCount \\\\ " >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "file size & $forsize & $revsize \\\\ " >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "\end{tabular}" >> ${mytex}.filestats
echo "\caption{\textbf{File stats}}" >> ${mytex}.filestats
echo "\end{table}" >> ${mytex}.filestats

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
    
    echo "$1"

    if [[ $1 == flu ]]; then
        echo "------> jhu Building Krona Graph..."
        kraken2krona.sh -i $sampleName-output.txt -k ${krakenDatabase} -o $sampleName-jhu-output.txt -r $sampleName-jhu-Krona_id_graphic.html
    else
        # Run Krona
        echo "------> ours Building Krona Graph..."
        cut -f2,3 $sampleName-output.txt > $sampleName-kronaInput.txt;
        /usr/local/bin/ktImportTaxonomy $sampleName-kronaInput.txt;
        mv taxonomy.krona.html $sampleName-Krona_identification_graphic.html;
        mv taxonomy.krona.html.files $sampleName-taxonomy.krona.html.files
    fi

    # Set variables and paths
    output=`ls *-output.txt`
    report=`ls *kraken_report.txt`
    ls -lh $output
    ls -lh $report
    
    printf "%s, %s file size, %s reads\n" ${forFile} ${forFileSize} ${forCount}

    if [[ $sampleType == "paired" ]]; then
        printf "%s, %s file size, %s reads\n" ${revFile} ${revFileSize} ${revCount}
    fi

    forCount=`grep -c '^+$' $forReads`
    echo "forCount: $forCount"
    declare -i x=${forCount}
    declare -i y=${revCount}
    echo "x: $x"
    echo "y: $y"
    echo "" | awk -v x=$x -v y=$y '{printf "Total single end read count: %'\''d\n", x+y}'

    #Section of results summary that calculates number of reads per type of organism (ex: ssRNA virus)
    echo "Summary of Kraken Findings"
    cRead=`grep -c "^C" $output`
    #uRead=`grep -c "^U" $output`
    virusreport=`awk ' $5 == "10239" {print $2}' $report`
    #let allReads=`wc -l $output | awk '{print $1}'`
    echo "allReads: $allReads"

    if [ -z $virusreport ]; then
        virusreport="zero"
    fi
    declare -i v=${virusreport}
    #declare -i z=${allReads}
    declare -i z=${forCount}
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
samtools sort ${refname}.raw.bam -o ${refname}.sorted.bam
echo "****Indexing Bam"
samtools index ${refname}.sorted.bam
rm ${refname}.sam
rm ${refname}.raw.bam
samtools view -h -b -F4 ${refname}.sorted.bam > ./$refname.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx2g -jar  ${picard} MarkDuplicates INPUT=${refname}.mappedReads.bam OUTPUT=$refname.dup.bam METRICS_FILE=$refname.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $refname.dup.bam"
samtools index $refname.dup.bam

echo "*** Calling VCFs with UnifiedGenotyper"
java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I $refname.dup.bam -o ${refname}.UG.vcf -nct 2

if [ -s ${refname}.UG.vcf ]; then
    echo "${g}.UG.vcf present continueing script"
else
    sleep 60
    echo "${refname}.UG.vcf is missing, try making again"
    java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${dupbam} -o ${refname}.UG.vcf -nct 2
    if [ -s ${refname}.UG.vcf ]; then
        echo "${refname}.UG.vcf present continueing script"
    else
        exit 1
    fi
fi

# make reference guided contig
java -Xmx2g -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${refname}.readreference.fasta -V ${refname}.UG.vcf

#echo ">${refname}" > ${refname}.readreference.temp; awk ' $8 ~ /^AN=2/ {print $4} ' ${refname}.UG.vcf | tr -d [:space:] >> ${refname}.readreference.temp; echo "" >> ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.readreference.fasta

echo "short BLAST"
blastn -query ${refname}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${refname}-readreference-max1-nt-id.txt -max_target_seqs 1 -outfmt "6 saccver"
#rm ${refname}.readreference.fasta

head -1 ${refname}-readreference-max1-nt-id.txt > ${refname}-readreference-max1-nt-id.txt.temp; mv ${refname}-readreference-max1-nt-id.txt.temp ${refname}-readreference-max1-nt-id.txt

if [ -s ${refname}-readreference-max1-nt-id.txt ]; then
    echo "Something was found"
else
    echo ">${refname}" > ${refname}.readreference.temp; awk ' $8 ~ /^AN=2/ {print $4} ' ${refname}.UG.vcf | tr -d [:space:] >> ${refname}.readreference.temp; echo "" >> ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.readreference.fasta
    blastn -query ${refname}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${refname}-readreference-max1-nt-id.txt -max_target_seqs 1 -outfmt "6 saccver"
    head -1 ${refname}-readreference-max1-nt-id.txt > ${refname}-readreference-max1-nt-id.txt.temp; mv ${refname}-readreference-max1-nt-id.txt.temp ${refname}-readreference-max1-nt-id.txt
    echo "No matches"
    #exit 1
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

python -u "${pythonGetFasta}" $acc $writeFile
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
        python -u "${pythonGetFasta}" $acc $writeFile
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        sleep 5
        if [ ! -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 120
            echo "Running python script to grab fasta file"
            python -u "${pythonGetFasta}" $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 320
                echo "Running python script to grab fasta file"
                python -u "${pythonGetFasta}" $acc $writeFile
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

echo "orgref: ${orgref}"
echo "refname: ${refname}"

grep '^#' ${orgref}-${refname}.hapreadyAll.vcf > header
grep -v '^#' ${orgref}-${refname}.hapreadyAll.vcf > snps

#chromname=`awk ' $1 !~ /^#/ {print $1}' ${orgref}-${refname}.hapreadyAll.vcf | head -1`
#if [ -s HCbody ]; then
#	echo "HCbody exists and has size greater than zero"
#	grep -v '^#' ${orgref}-${refname}.UG.vcf | awk -v c=$chromname 'BEGIN {FS="\t"; OFS="\t"} { if($10 == "./." ) print c, $2, $3, $4, "N", $6, $7, $8, "GT", "1"; else print c, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > UGbody
#else
#	echo "HCbody is empty"
#	grep -v '^#' ${orgref}-${refname}.UG.vcf | awk -v r=$refname 'BEGIN {FS="\t"; OFS="\t"} { if($10 == "./." ) print r, $2, $3, $4, "N", $6, $7, $8, "GT", "1"; else print r, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > UGbody
#fi

# The first input takes precedence
#cat HCbody UGbody | awk '{ if (a[$2]++ == 0) print $0; }' | sort -nk2,2 > body

# Get zero coverage regions the chromosome name from the UG positions needs to be updated to the HC names
awk -v c=$refname 'BEGIN{OFS="\t"}{if($10 == "./.") print c, $2, $3, $4, "N", $6, $7, $8, "GT", "1" }' ${orgref}-${refname}.UG.vcf > zeroformated

cat snps zeroformated | sort -nk2,2 > body

#cp ${orgref}-${refname}.hapreadyAll.vcf OLD-${orgref}-${refname}.hapreadyAll.vcf
cat header body > ${orgref}-${refname}.newhapreadyAll.vcf

# make reference guided contig THIS IS THE FINAL REFERNCE.  THE LAST ALIGNMENT WILL BE MADE AND BLAST RESULTS WILL BE BASED ON THIS CONSENSUS SEQUENCE
java -Xmx2g -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${orgref}-${refname}.reference_guided.fasta -V ${orgref}-${refname}.newhapreadyAll.vcf -IUPAC ${orgref}-${refname}


if [ $sampleType == "paired" ]; then
    java -Xmx2g -jar ${picard} SamToFastq INPUT=${orgref}-${refname}.dup.bam FASTQ=./${orgref}-${refname}-mapped_R1.fastq SECOND_END_FASTQ=./${orgref}-${refname}-mapped_R2.fastq
else
    java -Xmx2g -jar ${picard} SamToFastq INPUT=${orgref}-${refname}.dup.bam FASTQ=./${orgref}-${refname}-mapped_R1.fastq
fi

if [ $sampleType == "paired" ]; then
    mapCount=`cat ${orgref}-${refname}-mapped_R1.fastq ${orgref}-${refname}-mapped_R2.fastq | grep -c "^+$"`
else
    mapCount=`cat ${orgref}-${refname}-mapped_R1.fastq | grep -c "^+$"`
fi
echo "mapCount: $mapCount"

#Length of reference
countNTs=`grep -v ">" $ref | tr -d "\n" | wc | awk '{print $3}'`

#Number of nucleotides in reference with coverage
echo "*** Bamtools is getting coverage..."
bamtools coverage -in ${orgref}-${refname}.dup.bam | awk -v x=${orgref}-${refname} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${orgref}-${refname}-coveragefile

#Mean depth of coverage
meancov=`awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' ${orgref}-${refname}-coveragefile`

#Length of reference
#countNTs=`awk 'END{print NR}' ${orgref}-${refname}-coveragefile`

#count positions with coverage
covCount=`awk '{ if ($3 != 0) count++ } END { print count }' ${orgref}-${refname}-coveragefile`
echo "covCount $covCount"

declare -i x=${covCount}
declare -i y=${countNTs}

echo "----------------------> ${orgref} covCount: $x"
echo "----------------------> ${orgref} countNTs: $y"

#Percent of reference with coverage
perc=`awk -v x=$x -v y=$y 'BEGIN { print(x/y)*100}'`
echo "perc: $perc"

LC_NUMERIC=en_US
printf "%-45s %'10d %11.2f%% %'10dX\n" ${orgref}-${refname} $mapCount $perc $meancov >> ${summaryfile}.pre

echo "`printf "%'.0f\n" ${mapCount}` reads aligned to ${orgref}-${refname}"
echo "${perc}% genome coverage, $meancov"

awk -v x=${orgref}-${refname} 'BEGIN {OFS="\t"} {print x, $2, $3}' ${orgref}-${refname}-coveragefile | grep -v "segment_all" > ${orgref}-${refname}-samplecoveragefile

#rm *sorted.bam*
#rm *dup.bam*
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

#########################

java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${n}-${refname}.dict
    if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allow
if [[ $sampleType == "paired" ]]; then
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads $revReads > ${orgref}-${refname}.sam
else
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads > ${orgref}-${refname}.sam
fi

samtools view -bh -F4 -T $ref ${orgref}-${refname}.sam > ${orgref}-${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${orgref}-${refname}.raw.bam -o ${orgref}-${refname}.sorted.bam
echo "****Indexing Bam"
samtools index ${orgref}-${refname}.sorted.bam
rm ${orgref}-${refname}.sam
rm ${orgref}-${refname}.raw.bam
samtools view -h -b -F4 ${orgref}-${refname}.sorted.bam > ./${orgref}-${refname}.mappedReads.bam
samtools index ./${orgref}-${refname}.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx2g -jar  ${picard} MarkDuplicates INPUT=${orgref}-${refname}.mappedReads.bam OUTPUT=${orgref}-${refname}.dup.bam METRICS_FILE=${orgref}-${refname}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${orgref}-${refname}.dup.bam"
samtools index ${orgref}-${refname}.dup.bam

##############################
echo $ref

java -Xmx2g -jar ${gatk} -T ClipReads -R $ref -I ${orgref}-${refname}.dup.bam -o ${orgref}-${refname}.downsample.bam -filterNoBases -dcov 10
samtools index ${orgref}-${refname}.downsample.bam

rm *UG.vcf

echo "*** Calling VCFs with UnifiedGenotyper"
java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.UG.vcf -nct 2

if [ -s ${orgref}-${refname}.UG.vcf ]; then
    echo "${orgref}-${refname}.UG.vcf present continueing script"
else
    sleep 60
    echo "${orgref}-${refname}.UG.vcf is missing, try making again"
    java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.UG.vcf -nct 2
    if [ -s ${orgref}-${refname}.UG.vcf ]; then
        echo "${orgref}-${refname}.UG.vcf present continueing script"
    else
        exit 1
    fi
fi

if [[ $flu == yes ]]; then
        awk -v oref=$orgref 'BEGIN{OFS="\t"} $6 > 300 && $8 ~ /^AC=1/ {count++} END{print oref, count}' ${orgref}-${refname}.UG.vcf >> ${root}/ac1count
fi

#########

# make reference guided contig using Unified Genotyper
java -Xmx2g -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${refname}.readreference.fasta -V ${orgref}-${refname}.UG.vcf

echo ">${refname}" > ${refname}.readreference.temp; grep -v ">" ${refname}.readreference.fasta >> ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.fasta
rm ${refname}.readreference.*

# Align reads to Unified Genotyper made reference.  This is the last alignment that will make the final haplotypecaller reference guided assembly
ref="${refname}.fasta"
refname=`echo $ref | sed 's/\.fasta//'`
echo "ref: $ref and refname: $refname"

bwa index $ref
samtools faidx $ref
rm ${refname}.dict
java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=$ref OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${n}-${refname}.dict
    if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allow
if [[ $sampleType == "paired" ]]; then
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads $revReads > ${orgref}-${refname}.sam
else
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"${orgref}-${refname}""\t"PL:ILLUMINA"\t"PU:"${orgref}-${refname}"_RG1_UNIT1"\t"LB:"${orgref}-${refname}"_LIB1"\t"SM:"${orgref}-${refname}" $ref $forReads > ${orgref}-${refname}.sam
fi

samtools view -bh -F4 -T $ref ${orgref}-${refname}.sam > ${orgref}-${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${orgref}-${refname}.raw.bam -o ${orgref}-${refname}.sorted.bam
echo "****Indexing Bam"
samtools index ${orgref}-${refname}.sorted.bam
rm ${orgref}-${refname}.sam
rm ${orgref}-${refname}.raw.bam
samtools view -h -b -F4 ${orgref}-${refname}.sorted.bam > ./${orgref}-${refname}.mappedReads.bam
samtools index ./${orgref}-${refname}.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx2g -Xmx4g -jar  ${picard} MarkDuplicates INPUT=${orgref}-${refname}.mappedReads.bam OUTPUT=${orgref}-${refname}.dup.bam METRICS_FILE=${orgref}-${refname}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${orgref}-${refname}.dup.bam"
samtools index ${orgref}-${refname}.dup.bam
java -Xmx2g -jar ${gatk} -T ClipReads -R $ref -I ${orgref}-${refname}.dup.bam -o ${orgref}-${refname}.downsample.bam -filterNoBases -dcov 10
samtools index ${orgref}-${refname}.downsample.bam
########

#bam prepared now onto variant calling
java -Xmx2g -jar ${gatk} -R $ref -T HaplotypeCaller -ploidy 1 -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.hapreadyAll.vcf -bamout ${orgref}-${refname}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef 
java -Xmx2g -jar ${igvtools} index ${orgref}-${refname}.hapreadyAll.vcf

if [ -s ${orgref}-${refname}.hapreadyAll.vcf ]; then
    echo "${orgref}-${refname}.hapreadyAll.vcf present continueing script"
else
    sleep 60
    echo "${orgref}-${refname}.hapreadyAll.vcf is missing, try making again"
    java -Xmx2g -jar ${gatk} -R $ref -T HaplotypeCaller -ploidy 1 -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.hapreadyAll.vcf -bamout ${orgref}-${refname}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
    java -Xmx2g -jar ${igvtools} index ${orgref}-${refname}.hapreadyAll.vcf
    if [ -s ${orgref}-${refname}.hapreadyAll.vcf ]; then
        echo "${orgref}-${refname}.hapreadyAll.vcf present continueing script"
    else
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

python -u "${pythonGetFasta}" $acc $writeFile
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
        python -u "${pythonGetFasta}" $acc $writeFile
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        sleep 5
        if [ ! -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 120
            echo "Running python script to grab fasta file"
            python -u "${pythonGetFasta}" $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 320
                echo "Running python script to grab fasta file"
                python -u "${pythonGetFasta}" $acc $writeFile
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

#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow the most mismatch per read. -A [1] may be increased to increase the number of mismatches allowed
if [[ $sampleType == "paired" ]]; then
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"$refname""\t"PL:ILLUMINA"\t"PU:"$refname"_RG1_UNIT1"\t"LB:"$refname"_LIB1"\t"SM:"$refname" $ref $forReads $revReads > ${refname}.sam
else
bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"$refname""\t"PL:ILLUMINA"\t"PU:"$refname"_RG1_UNIT1"\t"LB:"$refname"_LIB1"\t"SM:"$refname" $ref $forReads > ${refname}.sam
fi

samtools view -bh -F4 -T $ref ${refname}.sam > ${refname}.raw.bam
echo "Sorting Bam"
samtools sort ${refname}.raw.bam -o ${refname}.sorted.bam
echo "****Indexing Bam"
samtools index ${refname}.sorted.bam
rm ${refname}.sam
rm ${refname}.raw.bam

echo "*** Calling VCFs with UnifiedGenotyper"
java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${refname}.sorted.bam -o ${refname}.UG.vcf -nct 2

if [ -s ${refname}.UG.vcf ]; then
    echo "${refname}.UG.vcf present continueing script"
else
    sleep 60
    echo "${refname}.UG.vcf is missing, try making again"
    java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${refname}.sorted.bam -o ${refname}.UG.vcf -nct 2
    if [ -s ${refname}.UG.vcf ]; then
        echo "${refname}.UG.vcf present continueing script"
    else
        exit 1
    fi
fi

# make reference guided contig
java -Xmx2g -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${refname}.readreference.fasta -V ${refname}.UG.vcf

sed 's/NN//g' < ${refname}.readreference.fasta > ${refname}.readreference.temp; mv ${refname}.readreference.temp ${refname}.readreference.fasta

echo "short BLAST"
blastn -query ${refname}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${refname}-readreference-max1-nt-id.txt -max_target_seqs 1 -outfmt "6 saccver"

rm *dict
rm *fasta*
rm *coveragefile
rm *vcf*
rm *bam*
rm ${refname}.readreference.fasta

}

#######################################################################################
#|||||||||||||||||||||| Function to make R graph of coverage ||||||||||||||||||||||||||
#######################################################################################
function plotR () {
cat > ./plotR.r << EOL
#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(scales)

arg <- commandArgs(trailingOnly=TRUE)

data <- read.csv(arg[1], header=FALSE, sep="\t")
names(data) <- c("Species", "position", "coverage")

pdf("myplot.pdf", width=20, height=6)

#ggplot(data, aes(x=position, y=log(coverage), colour=species, group=species)) + geom_point(size=2.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=10)))

#bp <- ggplot(data, aes(x=position, y=log10(coverage), colour=species, group=species)) + geom_point(size=2.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=10)))

#bp + scale_y_continuous(breaks=seq(0, 3.0, 0.5))

ggplot(data, aes(x=position, y=coverage, colour=Species, group=Species)) + geom_point(size=1.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=8)))+ scale_y_log10() + theme(axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"), axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"), title = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain")) + ylab("Depth of Coverage") + xlab("Position")

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

if [[ $flu == yes ]]; then

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
samtools sort ${refname}.raw.bam -o ${refname}.sorted.bam
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

    cp $targetref ./
    mkdir fastas
    cp $targetref ./fastas
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

	mv segment1 ./segment1_PB2/
        mv segment2 ./segment2_PB1/
        mv segment3 ./segment3_PA/
        mv segment4 ./segment4_HA/
        mv segment5 ./segment5_NP/
        mv segment6 ./segment6_NA/
        mv segment7 ./segment7_MP/
        mv segment8 ./segment8_NS/
	
	for s in segment*; do
		(cd $root
		cd $s
		echo "s: $s"
		pwd
		if [ "$bflag" ]; then
   		echo ""
    		echo " *** B FLAG ON, BUG FINDING MODE, SINGLE SAMPE PROCESSING *** "
    		echo ""
		for i in *fasta; do
                        cd ${root}/${s}
                        mkdir ${i%.fasta}
                        mv ${i} ${i%.fasta}
                        ln ${root}/*fastq* ${i%.fasta}
                        echo "working on $sampleName $s $i"
                        cd ${i%.fasta}; findbest
                done
		else
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
		fi
	wait
	cd ${root}/$s
	averagecov=`sed 's/%//g' ${root}/${s}/${sampleName}.findbest | sed 's/X$//' | sed 's/,//' | awk '{ sum += $2 } END { if (NR > 0) print sum / NR }'`
	meandepth=`sed 's/%//g' ${root}/${s}/${sampleName}.findbest | sed 's/X$//' | sed 's/,//' | awk '{ sum += $3 } END { if (NR > 0) print sum / NR }'`
	best=`sed 's/%//g' ${root}/${s}/${sampleName}.findbest | sed 's/X$//' | sed 's/,//' | awk -v x=$averagecov ' $2+1 > x {print $0}' | awk -v x=$meandepth ' $3+1 > x {print $0}' | sort -rk2,2 | head -1 | awk '{print $1}'`
	if [ -z $best  ]; then 
	
	best=`sed 's/%//' ${root}/${s}/${sampleName}.findbest | sed 's/X//' | sed 's/,//' | sort -nrk2,2 -nrk3,3 | head -1 | awk '{print $1}'` 
	
	fi
	echo "The best found: $best"
	echo "The best found: $best" >> ${root}/bestrefs.txt
	rm -r `ls | grep -v ${best}` 
	find . -name "*gz" -exec mv {} ./ \;
	find . -name "*fasta" -exec mv {} ./ \;
	segmentname=${PWD##*/}  # Segment name from working directory 
	rm -r ${best}
	mv *fasta ${segmentname}.fasta) &
                let count+=1
                [[ $((count%55)) -eq 0 ]] && wait
	done

	wait
	cd $root
	pwd

	for i in segment*; do 
		(echo ""; echo "####### $i ########"; echo ""; cd ${root}; cd $i; alignreads) &
                let count+=1
                [[ $((count%55)) -eq 0 ]] && wait
	done

else

if [ "$bflag" ]; then
    echo ""
    echo " *** B FLAG ON, BUG FINDING MODE, SINGLE SAMPE PROCESSING *** "
    echo ""
    for i in `ls $targetref`; do cp $i ./; done
    mkdir fastas
    cp $targetref ./fastas
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
    mkdir fastas
    cp $targetref ./fastas
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
printf "%-45s %10s %11s %10s\n" "reference used" "read count" "percent cov" "ave depth" >> ${summaryfile}
sort -k1,1 < ${summaryfile}.pre >> ${summaryfile}

echo "" >> ${mytex}.alignmentstats
echo "\vspace{5mm}" >> ${mytex}.alignmentstats
echo "" >> ${mytex}.alignmentstats

echo "\begin{table}[H]" >> ${mytex}.alignmentstats
echo "\begin{tabular}{ l | l | l | l }" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
echo "reference used & read count & percent cov & ave depth \\\\" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4}' ${summaryfile}.pre | sort -k1,1 | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\ \\hline:' | sed 's/[_]/\\_/g' | sed 's/[%]/\\%/g' >> ${mytex}.alignmentstats
echo "\end{tabular}" >> ${mytex}.alignmentstats
echo "\caption{\textbf{Alignment stats (reference guided)}}" >> ${mytex}.alignmentstats
echo "\end{table}" >> ${mytex}.alignmentstats

currentdate=`date`
sort -k1,1 < ${summaryfile}.pre | awk -v name="$sampleName" -v curdate="$currentdate" 'BEGIN{OFS="\t"} {print curdate, name, $1, $2, $3, $4}' >> $idscriptrunsummary

sort -k1,1 < ${summaryfile}.pre | awk -v name="$sampleName" -v curdate="$currentdate" 'BEGIN{OFS="\t"} {print curdate, name, $1, $2, $3, $4}' >> /bioinfo11/TStuber/Results/viruses/idvirus_run_summary.txt

echo "" >> ${emailbody}
echo "Alignment stats (reference guided):" >> ${emailbody}
printf "%-45s %10s %11s %10s\n" "reference used" "read count" "percent cov" "ave depth" >> ${emailbody}
sort -k1,1 < ${summaryfile}.pre >> ${emailbody}

#####################

for i in `find . -name "*samplecoveragefile"`; do
    cat $i >> allsamplecoveragefile
done
pwd

plotR allsamplecoveragefile $sampleName
mv myplot.pdf ${sampleName}.assembly_graph.pdf
pwd

cd $root

pwd 

mkdir ${sampleName}-reference_guided_assemblies
for i in `find . -name "*reference_guided.fasta"`; do
    cp $i ${sampleName}-reference_guided_assemblies
done
pwd

cd ${root}/${sampleName}-reference_guided_assemblies
pwd

for i in `ls *fasta | sort -k1,1`; do
    echo ">${i%.reference_guided.fasta}" >> ${sampleName}.consensus.reads.fasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | awk 'length > x { x = length; y = $0 } END { print y }' >> ${sampleName}.consensus.reads.fasta
done
pwd

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
java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${refname}.dict

if [ -s ${ref}.fai ] && [ -s ${refname}.dict ]; then
echo "Index and dict are present, continue script"
else
sleep 5
echo "Either index or dict for reference is missing, try making again"
samtools faidx $ref
java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${n}-${refname}.dict
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
samtools sort ${orgref}-${refname}.raw.bam -o ${orgref}-${refname}.sorted.bam
echo "****Indexing Bam"
samtools index ${orgref}-${refname}.sorted.bam
rm ${orgref}-${refname}.sam
rm ${orgref}-${refname}.raw.bam
samtools view -h -b -F4 ${orgref}-${refname}.sorted.bam > ./${orgref}-${refname}.mappedReads.bam
samtools index ./${orgref}-${refname}.mappedReads.bam

echo "***Marking Duplicates"
java -Xmx2g -jar  ${picard} MarkDuplicates INPUT=${orgref}-${refname}.mappedReads.bam OUTPUT=${orgref}-${refname}.dup.bam METRICS_FILE=${orgref}-${refname}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${orgref}-${refname}.dup.bam"
samtools index ${orgref}-${refname}.dup.bam
java -Xmx2g -jar ${gatk} -T ClipReads -R $ref -I ${orgref}-${refname}.dup.bam -o ${orgref}-${refname}.downsample.bam -filterNoBases -dcov 10
samtools index ${orgref}-${refname}.downsample.bam

# Make a quick and simple VCF to highlight possible problem areas of the consensus
java -Xmx2g -jar ${gatk} -R $ref -T UnifiedGenotyper -glm BOTH -I ${orgref}-${refname}.downsample.bam -o ${orgref}-${refname}.UG.vcf -nct 8

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

wait
# Go back to the assemble folder
cd $assemblyfolder
wait
###########

contigcount=`grep -c ">" ${sampleName}.consensus.reads.fasta`
sed 's/NNN//g' ${sampleName}.consensus.reads.fasta > ${sampleName}.consensusnoN.reads.fasta
pwd

echo "nt BLAST $contigcount contigs..."
blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${sampleName}-consensus-max1-nt.txt -max_target_seqs 1 -outfmt "6 qseqid qlen slen pident mismatch evalue bitscore stitle saccver"


echo "" >> ${summaryfile}

hsegment=`grep "segment4" ${sampleName}-consensus-max1-nt.txt | sed 's/.*\(H[0-9]\{1,2\}\).*/\1/' | head -1`
echo "hsegment: $hsegment"
nsegment=`grep "segment6" ${sampleName}-consensus-max1-nt.txt | sed 's/.*\(N[0-9]\{1,2\}\).*/\1/' | head -1`
echo "nsegment: $nsegment"

#If matching segment is labeled "mixed" H and N subtypes are not found above and the entire line identification is found
#therefore if characters in variable are > 3 variable is changed to "mixed"
if [[ ${#hsegment} -gt 3 ]]; then
        echo "hsegment is greater thatn 3, likely mixed.  Changing variable to mixed"
        hsegment="H?"
fi

if [[ ${#nsegment} -gt 3 ]]; then 
	echo "nsegment is greater thatn 3, likely mixed.  Changing variable to mixed"
	nsegment="N?"
fi

subtype=${hsegment}${nsegment}

# echo subtype to terminal
i=0; while [ $i -lt 11 ]; do echo "*** Subtype: $subtype"; i=$[$i+1]; done

echo "subtype: $subtype"

if [[ -n $subtype ]]; then
        echo "Subtype: $subtype" >> ${summaryfile}
	echo "Subtype: $subtype" >> $idscriptrunsummary
	echo "Subtype: $subtype" >> /bioinfo11/TStuber/Results/viruses/idvirus_run_summary.txt

	sed -i "s/XXXXXHNTYPEXXXXXXX/Subtype: $subtype/" $mytex

else
	sed -i "s/XXXXXHNTYPEXXXXXXX/$argUsed/" $mytex

fi

echo "--------------------------------------------------" >> ${summaryfile}
echo "*** NT database ***" >> ${summaryfile}


echo "" >> $mytex
echo "\vspace{5mm}" >> $mytex
echo "" >> $mytex

echo "\begin{table}[H]" >> $mytex
echo "\tiny" >> $mytex
echo "\begin{tabular}{ p{4cm} | l | l | l | l | l | l | p{4cm} }" >> $mytex
echo "\hline" >> $mytex
echo "ID & qlength & slength & \\% id & mis & evalue & bscore & Description \\\\" >> $mytex
echo "\hline" >> $mytex
echo "\hline" >> $mytex
cut -f1-8 ${sampleName}-consensus-max1-nt.txt | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\ \\hline:' | sed 's/_/\\_/g' >> $mytex
echo "\end{tabular}" >> $mytex
echo "\caption{\textbf{NT database}}" >> $mytex
echo "\end{table}" >> $mytex

awk 'BEGIN{printf "%-45s %-8s %-8s %-8s %-3s %-6s %-6s %-1s\n", "query ID", "qlength", "slength", "% id", "mis", "evalue", "bscore", "Description"} {printf "%-45s %-8s %-8s %-8s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' ${sampleName}-consensus-max1-nt.txt >> ${summaryfile}
echo "" >> ${summaryfile}

sort -k1,1 < ${summaryfile}.pre > ${summaryfile}.sorted
echo "*** ${sampleName} ${subtype}" >> /scratch/report/idemailsummary
paste ${summaryfile}.sorted ${sampleName}-consensus-max1-nt.txt | awk 'BEGIN{printf "%-41s %-11s %-8s %-10s %-3s %-6s %-6s %-1s\n", "ID", "read count", "per cov", "ave depth", "mis", "evalue", "bscore", "Description"} {printf "%-41s %-11s %-8s %-10s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' >> /scratch/report/idemailsummary

#awk 'BEGIN{printf "%-41s %-11s %-8s %-10s %-3s %-6s %-6s %-1s\n", "ID", "read count", "per cov", "ave depth", "mis", "evalue", "bscore", "Description"} {printf "%-41s %-11s %-8s %-10s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' ${summaryfile}.sorted  >> /scratch/report/idemailsummary

rm ${summaryfile}.pre ${summaryfile}.sorted

pingyrdb=`echo ${subtype} | egrep -m 1 -o "H5N1|H5N2|H5N8"`
echo "In the pingyrdb varable: $pingyrdb"

if [[ -n $pingyrdb ]]; then
    echo "pintail-gyrfalcon BLAST $contigcount contigs..."
    blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/blastdb-pintail-gyrfalcon/pintail-gyrfalcon.fsa -num_threads 20 -out ${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt -outfmt 1
    blastn -query ${sampleName}.consensusnoN.reads.fasta -db /data/BLAST/blastdb-pintail-gyrfalcon/pintail-gyrfalcon.fsa -num_threads 20 -out ${sampleName}-consensus-fmt6-pintail-gyrfalcon.txt -outfmt "6 qseqid qlen slen pident mismatch evalue bitscore stitle saccver"
    echo "--------------------------------------------------" >> ${summaryfile}
    echo "*** Pintail/Gyrfalcon database ***" >> ${summaryfile}

    awk 'BEGIN{printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s\n", "query ID", "qlength", "slength", "% id", "mis", "evalue", "bscore", "Description"} {printf "%-30s %-8s %-8s %-8s %-3s %-6s %-6s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s %-1s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27}' ${sampleName}-consensus-fmt6-pintail-gyrfalcon.txt >> ${summaryfile}

echo "" >> $mytex
echo "\vspace{5mm}" >> $mytex
echo "" >> $mytex

echo "\begin{table}[H]" >> $mytex
echo "\tiny" >> $mytex
echo "\begin{tabular}{ p{4cm} | l | l | l | l | l | l | p{4cm} }" >> $mytex
echo "\hline" >> $mytex
echo "ID & qlength & slength & \\% id & mis & evalue & bscore & Description \\\\" >> $mytex
echo "\hline" >> $mytex
echo "\hline" >> $mytex
cut -f1-8 ${sampleName}-consensus-fmt6-pintail-gyrfalcon.txt | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\ \\hline:' | sed 's/_/\\_/g' >> $mytex
echo "\end{tabular}" >> $mytex
echo "\caption{\textbf{Pintail/Gyrfalcon database}}" >> $mytex
echo "\end{table}" >> $mytex

fi

###########################

mkdir ${root}/share_folder
cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/share_folder

############################
cd ${root}/share_folder

#Format fasta headers for NCBI

if [[ ${genotypingcodes} == "NEED TO SET" ]]; then
    echo "genotyping codes not given"
    cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/${sampleName}-submissionfile.fasta
    sed -i "s/XXXXXSTRAINNAMEXXXXXXX//" $mytex	
else
    echo "sampleName: $sampleName"

    grep "$sampleName" $genotypingcodes | head -1 > ${sampleName}.information
    if [[ -s ${sampleName}.information ]]; then
        echo "file exists and is greater than zero, continue"
        #column 1: sample
        sample=`awk 'BEGIN{FS="\t"}{print $1}' ${sampleName}.information`
        echo "sample $sample"

        #column 2: species
        species=`awk 'BEGIN{FS="\t"}{print $2}' ${sampleName}.information`
        speciesspace=`awk 'BEGIN{FS="\t"}{print $2}' ${sampleName}.information | sed 's/_/ /g'`
    	echo "species $species"
        echo "speciesspace $speciesspace"

        #column 3: state
        state=`awk 'BEGIN{FS="\t"}{print $3}' ${sampleName}.information`
        echo "state $state"
        statespace=`awk 'BEGIN{FS="\t"}{print $3}' ${sampleName}.information | sed 's/_/ /g'`
        syear=`echo "$sample" | sed 's/-.*//'`
        echo "syear $syear"
        sampleyear=`echo "20${syear}"`
        echo "sampleyear $sampleyear"

        sed 's/>.*seg.*1_.*/>Seq1/' ${sampleName}.consensus.reads.fasta | sed 's/>.*seg.*2_.*/>Seq2/' | sed 's/>.*seg.*3_.*/>Seq3/' | sed 's/>.*seg.*4_.*/>Seq4/' | sed 's/>.*seg.*5_.*/>Seq5/' | sed 's/>.*seg.*6_.*/>Seq6/' | sed 's/>.*seg.*7_.*/>Seq7/' | sed 's/>.*seg.*8_.*/>Seq8/' > ${sampleName}.temp

#noyear=`echo $sample | sed -e "s/$syear-//"`
#echo "noyear $noyear"

	sed -i "s:XXXXXSTRAINNAMEXXXXXXX:A/${speciesspace}/${statespace}/${sample}/${sampleyear}:" $mytex
        	
	# Create "here-document"
cat >./param.txt <<EOL

>Seq1 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 1, polymerase PB2 (PB2) gene, complete cds.
>Seq2 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 2, polymerase PB1 (PB1) gene, complete cds.
>Seq3 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 3, polymerase PA (PA) gene, complete cds.
>Seq4 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 4, hemagglutinin (HA) gene, complete cds.
>Seq5 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 5, nucleoprotein (NP) gene, complete cds.
>Seq6 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 6, neuraminidase (NA) gene, complete cds.
>Seq7 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 7, matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds.
>Seq8 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 8, non-structural protein NS1 and non-structural protein NS2 (NS) gene, complete cds.

EOL

        awk 'NR==FNR{a[$1]=$0;next} ($1 in a){ print a[$1]; next}1' param.txt ${sampleName}.temp > ${root}/${sampleName}-submissionfile.fasta

    else
	sed -i "s/XXXXXSTRAINNAMEXXXXXXX//" $mytex
	sed -i "s/XXXXXHNTYPEXXXXXXX/$argUsed/" $mytex
        echo "metadata not available"
        cp ${root}/${sampleName}-reference_guided_assemblies/${sampleName}.consensus.reads.fasta ${root}/${sampleName}-submissionfile.fasta
    fi
pwd

#rm *temp
rm *information
rm param.txt

fi

###########################
echo "Making IRD for: $sampleName"

# Create "here-document"
cat >./ird_param.txt <<EOL

>Seq1 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.PB2
>Seq2 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.PB1
>Seq3 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.PA
>Seq4 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.HA
>Seq5 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.NP
>Seq6 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.NA
>Seq7 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.MP
>Seq8 Unique_Sample_Identifier:${sampleName}|Unique_Sequence_Identifier:${sampleName}.NS

EOL

awk 'NR==FNR{a[$1]=$0;next} ($1 in a){ print a[$1]; next}1' ird_param.txt ${sampleName}.temp | sed 's/>Seq../>/' > ${root}/${sampleName}-IRD-submissionfile.fasta
pwd

rm *temp
rm *information
rm param.txt

###########################

cd $root
cp ${sampleName}.assembly_graph.pdf graphic.pdf

echo "\vspace{5mm}" >> $mytex
echo "" >> $mytex

echo "\begin{figure}[H]" >> $mytex
echo "\begin{flushleft}" >> $mytex
echo "\textbf{Coverage Graph}\par\medskip" >> $mytex 
echo "\end{flushleft}" >> $mytex

echo "\includegraphics[width=450pt]{graphic.pdf}" >> $mytex
echo "" >> $mytex
echo "\end{figure}" >> $mytex
echo "" >> $mytex

#add file and alignment stats
cat ${mytex}.filestats >> $mytex
cat ${mytex}.alignmentstats >> $mytex

#---
# if flu then add a table to report that shows the if suspected mix segment
if [[ $flu == yes ]]; then
        awk 'BEGIN {OFS="\t"} {if ($2 > 10) print $1, "mix"; else print $1, "pure"}' ${root}/ac1count | sort > ${root}/sortedac1count
        echo "" >> ${mytex}.acstats
        echo "\vspace{5mm}" >> ${mytex}.acstats
        echo "" >> ${mytex}.acstats

        echo "\begin{table}[H]" >> ${mytex}.acstats
        echo "\begin{tabular}{ l | l }" >> ${mytex}.acstats
        echo "\hline" >> ${mytex}.acstats
        echo "segment & findings \\\\" >> ${mytex}.acstats
        echo "\hline" >> ${mytex}.acstats 
        echo "\hline" >> ${mytex}.acstats
        awk 'BEGIN{OFS="\t"}{print $1, $2}' ${root}/sortedac1count | sort -k1,1 | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\ \\hline:' | sed 's/[_]/\\_/g' | sed 's/[%]/\\%/g' >> ${mytex}.acstats
        echo "\end{tabular}" >> ${mytex}.acstats
        echo "\caption{\textbf{Mix Test Results}}" >> ${mytex}.acstats
        echo "\end{table}" >> ${mytex}.acstats

        cat ${mytex}.acstats >> $mytex
fi
#---

echo "\end{document}" >> $mytex
echo "" >> $mytex

pdflatex $mytex
pdflatex $mytex
mv $sampleName.pdf ${sampleName}-report.pdf

#rm *fastq*
echo "" >> ${summaryfile}
echo "" >> ${summaryfile}
#echo "Files copied to: $bioinfoVCF" >> ${summaryfile}
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
cp ${sampleName}-reference_guided_assemblies/${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt ${root}

#enscript ${summaryfile} -B -j -r -f "Courier5" -o - | ps2pdf - ${sampleName}-report.pdf

if [ -e ${root}/${sampleName}-report.pdf ]; then
	ls ${root}/${sampleName}-report.pdf > emailfiles
fi

if [ -e ${root}/${sampleName}-submissionfile.fasta ]; then
	ls ${root}/${sampleName}-submissionfile.fasta >> emailfiles
fi

if [ -e ${root}/${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt ]; then
	ls ${root}/${sampleName}-consensus-blast_alignment-pintail-gyrfalcon.txt >> emailfiles
fi

if [ -e ${root}/$sampleName-jhu-Krona_id_graphic.html ]; then 
	ls ${root}/$sampleName-jhu-Krona_id_graphic.html >> emailfiles
fi

if [ -e ${root}/$sampleName-Krona_identification_graphic.html ]; then 
	ls ${root}/$sampleName-Krona_identification_graphic.html >> emailfiles
fi

if [ -e ${root}/kraken/${sampleName}-kraken_report.txt ]; then
        ls ${root}/kraken/${sampleName}-kraken_report.txt >> emailfiles
fi

if [[ $sampleType == "paired" ]]; then
	echo "paried data, not checking for C insert"
else
	if [[ -n $pingyrdb ]]; then
		noc=`egrep -c "GAGTTGACATAAACCAGGCCACGC|GCGTGGCCTGGTTTATGTCAACTC" $forReads`
		cinsert=`egrep -c "GAGTTGACATAAACCCAGGCCACGC|GCGTGGCCTGGGTTTATGTCAACTC" $forReads`
        insert1=`egrep -c "GAGTTGACATAAA[AGT]CCAGGCCACGC|GCGTGGCCTGG[ACT]TTTATGTCAACTC" $forReads`
        	insertA=`egrep -c "GAGTTGACATAAAACCAGGCCACGC|GCGTGGCCTGGTTTTATGTCAACTC" $forReads`
		insertT=`egrep -c "GAGTTGACATAAATCCAGGCCACGC|GCGTGGCCTGGATTTATGTCAACTC" $forReads`
		insertG=`egrep -c "GAGTTGACATAAAGCCAGGCCACGC|GCGTGGCCTGGCTTTATGTCAACTC" $forReads`

	insert2=`egrep -c "GAGTTGACATAAAC[AGT]CAGGCCACGC|GCGTGGCCTG[ACT]GTTTATGTCAACTC" $forReads`
        insert3=`egrep -c "GAGTTGACATAAACC[AGT]AGGCCACGC|GCGTGGCCT[ACT]GGTTTATGTCAACTC" $forReads`
		echo "" >> ${emailbody}
		echo "POSITION 446 read counts:" >> ${emailbody}
		echo "(-CC) No insert count: $noc" >> ${emailbody}
		echo "(CCC) C insert count: $cinsert" >> ${emailbody}
        echo "" >> ${emailbody}
	echo "([AGT]CC) A,G,T insert count at position 1: $insert1" >> ${emailbody}
        	echo "           (ACC) A insert count at position 1: $insertA" >> ${emailbody}
		echo "           (TCC) T insert count at position 1: $insertT" >> ${emailbody}
		echo "           (GCC) G insert count at position 1: $insertG" >> ${emailbody}	
	echo "" >> ${emailbody}
	echo "(C[AGT]C) A,G,T insert count at position 2: $insert2" >> ${emailbody}
        echo "(CC[AGT]) A,G,T insert count at position 3: $insert3" >> ${emailbody}
		echo "" >> ${emailbody}
	else
		echo "pingyrdb not being referenced, therefore not checking for C insert"
	fi
fi

rm *fastq*

rm ${sampleName}.assembly_graph.pdf
rm ${sampleName}.aux
rm ${sampleName}.log
rm ${sampleName}.tex.alignmentstats
rm ${sampleName}.tex.filestats
rm allsamplecoveragefile
rm bestrefs.txt
rm writelist

#Cleanup
rm -r `ls | egrep -v "$myfile|${myfile.tex}.pdf|kraken|emailfile|emailfiles|bestrefs.txt|$0|igv_alignment|originalreads|summaryfile|report.pdf|_graphic.html|-consensus-blast_alignment-pintail-gyrfalcon.txt|-submissionfile.fasta|assembly_graph.pdf"`

pwd > ./fastas/filelocation.txt

echo "$sampleName"
echo "$subtype"
echo "$argUsed"

if [ "$eflag" ]; then
	# eflag is used when script is called from idemail.sh
	# making summary file to send in email
	echo "Files copied to: ${bioinfoVCF}" >> /scratch/report/idemailsummary
	echo "" >> /scratch/report/idemailsummary
	rm emailfile*
	echo "Copying to ${bioinfoVCF}"
        cp -r $PWD ${bioinfoVCF} &
else
	# else when idvirus.sh is ran on its own
	if [ "$mflag" ]; then
    		email_list="tod.p.stuber@usda.gov"
    		cat ${emailbody} | mutt -s "Sample: ${sampleName}, $subtype Reference_Set: $argUsed" -a `cat emailfiles` -- $email_list
    		rm emailfiles
	else
    		echo "" >> ${emailbody}
    		echo "Files copied to: ${bioinfoVCF}" >> ${emailbody}		
    		cat ${emailbody} | mutt -s "Sample: ${sampleName}, $subtype Reference_Set: $argUsed" -a `cat emailfiles` -- $email_list

    		rm ${emailbody}
    		rm emailfiles
    		echo "Copying to ${bioinfoVCF}"
    		cp -r $PWD ${bioinfoVCF} &
	fi
fi

rm
echo "****************************** END ******************************"
pwd


# Created: 2015-06-12, tstuber
