#!/bin/sh

: <<'END'
This script is the second script in a two script workflow.  Script 2 genotypes Mycobacterium tuberculosis complex and Brucella species from SNP data contained in VCFs.  It operates on VCFs generated with the same reference output from script 1.  VCFs are collected into a single working directory.  Comparisons are output as SNP tables and alignment files to view as trees in your program of choice.

Script 2 will run and output tables and alignments from just the data generated from script 1, however the data will be more informative if additional files are provide.  Those files are:
1) A file that contains positions to cluster individual isolates/VCFs into groups, subgroups and clades.
2) Files that contain positions to remove from the analysis.

Paradigm
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show parsimony informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

END
echo ""
echo "****************************** START ******************************"
echo ""

echo "Start Time: `date`" > sectiontime
starttime=`date +%s`
argUsed="$1"
uniqdate=`date "+%Y-%m-%dat%Hh%Mm%Ss"`
dircalled=`pwd`
echo "start time: $uniqdate"

# Set flags
# flag -c with look for positions to filter.  By default, with no -c, this will not be done.
# flag -m will email just "M"e
# flag -e will run the bovis "E"lite representative samples
# flag -a get "a"ll_vcf alignment table

cflag=
mflag=
eflag=
aflag=
while getopts 'cmea' OPTION; do
    case $OPTION in
        c) cflag=1
        ;;
        m) mflag=1
        ;;
        e) eflag=1
        ;;
        a) aflag=1
        ;;
        ?) echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))

####################################################
filterdir="/home/shared/${uniqdate}-FilterFiles"
mkdir ${filterdir}
FilterDirectory=${filterdir} #Files containing positions to filter
####################################################

function filterFileCreations () {

# Use to make filter files from the text pasted from the Excel worksheet.
# working directory does not need to be set.
#   Set variables:

# Path to txt file containing paste from Excel worksheet.
filterFile="${filterdir}/filterFile.txt"

# Number of columns in Excel worksheet
columns=`head $filterFile | awk 'BEGIN{ FS="\t"; OFS="\t" }  END {print NF}'`

# Location filter files are output to.
output="${filterdir}"

let columns=columns+1
rm ${output}*
echo "Filter sets: $columns"
echo "Extracting from Excel to text files..."

count=1
while [ $count -lt ${columns} ]; do
    #echo ${count}
    filename=`awk -v x=$count 'BEGIN{FS=OFS="\t"}{print $x}' $filterFile | head -n1`
    #echo "Filename: $filename"
    awk -v x=$count 'BEGIN{FS=OFS="\t"} FNR>1 {print $x}' $filterFile | grep -v "^$" > ${output}/${filename}.list
    let count=count+1
done
rm $filterFile
for i in ${output}/*.list; do
    (base=`basename "$i"`
    readyfile=`echo $base | sed 's/\..*//'`

    touch ${output}/${readyfile}.txt

    mylist=`cat $i`

    for l in $mylist; do
        pos1=`echo $l | sed 's/-/ /g' | awk '{print $1}'`
        pos2=`echo $l | sed 's/-/ /g' | awk '{print $2}'`
        #echo $pos2
            if [[ -z "$pos2" ]]
            then
            let pos2=pos1+1
                while [ $pos1 -lt $pos2 ]; do
                #echo $pos1
                echo $pos1 >> ${output}/${readyfile}.txt
                let pos1=pos1+1
                done
            else
            let pos2=pos2+1
                while [ $pos1 -lt $pos2 ]; do
                #echo $pos1
                echo $pos1 >> ${output}/${readyfile}.txt
                let pos1=pos1+1
                done
            fi
        done) &
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

rm ${output}/*.list
}
####################################################
function parseXLS () {
# Create "here-document"
cat >./inputXLS.py <<EOL
#!/usr/bin/env python

import os
import xlrd
from sys import argv

script, input = argv

wb = xlrd.open_workbook(input)
wb.sheet_names()
#sh = wb.sheet_by_index(1)
sh = wb.sheet_by_name(u'New groupings')
for rownum in range(sh.nrows):
    print sh.row_values(rownum)

EOL

chmod 755 ./inputXLS.py

./inputXLS.py $excelinfile

rm ./inputXLS.py

}
#####################################################

# Environment controls:

if [[ $1 == ab1 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/abortus1/script_dependents/Abortus1_Defining_SNPs.txt"
    #coverageFiles="/bioinfo11/TStuber/Results/brucella/Abortus1/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/abortus1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/abortus1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/abortus1/vcfs"
    echo "vcftofasta.sh ran as Brucella abortus bv 1, 2 or 4"
    echo "Script vcftofasta.sh ran using Brucella abortus bv 1, 2 or 4 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov, Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == mel ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/melitensis/script_dependents/Mel_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/melitensis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/melitensis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/melitensis/vcfs"
    echo "vcftofasta.sh ran as B. melitensis"
    echo "Script vcftofasta.sh ran using B. melitensis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == suis1 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/suis1/script_dependents/Suis1_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/suis1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/suis1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/suis1/vcfs"
    echo "vcftofasta.sh ran as B. suis bv1"
    echo "Script vcftofasta.sh ran using B. suis bv1 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == suis2 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/suis2/script_dependents/suis2_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/suis2/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/suis2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/suis2/vcfs/"
    echo "vcftofasta.sh ran as B. suis bv2"
    echo "Script vcftofasta.sh ran using B. suis bv2 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == suis3 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/suis3/script_dependents/Suis3_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/suis3/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/suis3/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/suis3/vcfs"
    echo "vcftofasta.sh ran as B. suis bv3"
    echo "Script vcftofasta.sh ran using B. suis bv3 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == suis4 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/suis4/script_dependents/Suis4_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/suis4/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/_Brucela/suis4/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/suis4/vcfs"
    echo "vcftofasta.sh ran as B. suis bv4"
    echo "Script vcftofasta.sh ran using B. suis bv4 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == canis ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/canis/script_dependents/Canis_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/canis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/canis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/canis/vcfs"
    echo "vcftofasta.sh ran as B. canis"
    echo "Script vcftofasta.sh ran using B. canis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"


elif [[ $1 == ceti1 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/ceti1/script_dependents/Ceti1_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/ceti1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/ceti1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/ceti1/vcfs"
    echo "vcftofasta.sh ran as B ceti group 1"
    echo "Script vcftofasta.sh ran using B ceti group 1 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"


elif [[ $1 == ceti2 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/ceti2/script_dependents/Ceti2_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/ceti2/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/ceti2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/ceti2/vcfs"
    echo "vcftofasta.sh ran as B ceti group 2"
    echo "Script vcftofasta.sh ran using B ceti group 2 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"


elif [[ $1 == ovis ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/ovis/script_dependents/Ovis_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/bioinfo11/TStuber/Results/brucella/ovis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/ovis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/ovis/vcfs"
    echo "vcftofasta.sh ran as B. ovis"
    echo "Script vcftofasta.sh ran using B. ovis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

elif [[ $1 == bovis ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/comparisons"
    echo "vcftofasta.sh ran as M. bovis"
    echo "Script vcftofasta.sh ran using M. bovis variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    if [ "$eflag" ]; then
        echo "Only the "elite" bovis isolates are being ran"
    else
        echo "All bovis are being ran"
        echo "Like to run selected isolates? Use... vcftofasta.sh -e bovis"
    fi

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/Filtered_Regions.xlsx
    # Excel tab label "New groupings"

    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb1 ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb1/tb1DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb1/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb1/tb1Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  >${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb2 ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb2/tb2DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb2/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb2/tb2Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb3 ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb3/tb3DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb3/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb3/tb3Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb4a ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4a/tb4aDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4a/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4a/tb4aFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb4b ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4b/tb4bDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4b/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4b/tb4bFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb5 ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb5/tb5DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb5/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb5/tb5Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == tb6 ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/Untitled.tab"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb6/tb6DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb6/vcfs"
    echo "vcftofasta.sh ran as ${1}"
    echo "Script vcftofasta.sh ran using ${1} variables" >> section5
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb6/tb6Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == para ]]; then
    genotypingcodes="/bioinfo11/TStuber/Results/mycobacterium/mac/tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/vcfs"
    echo "vcftofasta.sh ran as M. paraTB"
    echo "Script vcftofasta.sh ran using para variables" >> section5
    email_list="tod.p.stuber@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov Boojala.Vijay.Reddy@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"

    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/vcfs/Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

elif [[ $1 == h5n2 ]]; then

	genotypingcodes="/bioinfo11/MKillian/Analysis/results/snp-genotypingcodes.txt"
	# This file tells the script how to cluster VCFs
	DefiningSNPs="/bioinfo11/MKillian/Analysis/results/influenza/h5n2/snp_analysis/script2/Defining_SNPs_H5N2.txt"
	FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
	FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
	FilterDirectory="/bioinfo11/MKillian/Analysis/results/influenza/h5n2/snp_analysis/script2/FilterFiles" #Files containing positions to filter
	RemoveFromAnalysis="bioinfo11/TStuber/Results/mycobacterium/vcfs/RemoveFromAnalysis.txt"
	QUAL=300 # Minimum quality for calling a SNP
	lowEnd=1
	highEnd=350 # QUAL range to change ALT to N
	bioinfoVCF="/bioinfo11/MKillian/Analysis/results/influenza/h5n2/snp_analysis/script2/"
	echo "vcftofasta.sh ran as H5N2"
	echo "Script vcftofasta.sh ran using h5n2 variables" > section5
	email_list="tod.p.stuber@usda.gov" #Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov
	#for i in *vcf; do awk 'BEGIN{OFS="\t"}$1 ~ /seg1/ || $1 ~ /^#/ {print $0}' $i > ../h5n2_2015-10-03-seg1/${i%.vcf}-seg1.vcf; done

else

    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, bovis, tb1, tb2, tb3, tb4a, tb4b, tb5, tb6, para, h5n2"
    echo ""
    echo "Set optional flags"
    echo "flag -c with look for positions to filter.  By default, with no -c, this will not be done."
    echo "flag -m will email just "M"e"
    echo "flag -e will run the bovis "E"lite representative samples"
    echo "flag -a get "a"ll_vcf alignment table"
    echo ""
    echo "Example: [prompt]$ vcftofasta.sh -mea bovis"
    echo ""
    rm sectiontime
    exit 1

fi
#################################################################################
# Set variables:

# Sed searches put into variables
tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9,FM]\{4,6\}\).*/\1/' #Only tb Number, *laboratory specific*
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file

NR_CPUS=50 # Computer cores to use when analyzing
LIMIT_CPUS=10

Ncov=1 # Coverage below this value will be changed to -

fulDir=$PWD # Current working directory, do not change.

# Remove selected files from comparison
# Use file:  /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysis.txt

function removeIsolates () {

cat ${RemoveFromAnalysis} | tr '\r' '\n' | awk '{print $1}' > /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysisUnixReady.txt

removeList=`cat /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysisUnixReady.txt`

for i in $removeList; do
    rm *${i}* > /dev/null 2>&1
done

rm /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysisUnixReady.txt

}

#################################################################################

# If there are 2 vcf files with the same name one of the files might unknowingly
# get cut out of the analysis and keep the undesired vcf instead.  This will
# alert if 2 vcf with the same TB number are present.
# The regular expression used in sed should be changed based on vcf naming convention

function testDuplicates () {

echo "Checking for empty or duplicate VCFs."

directorytest="${PWD##*/}"
	if [[ $directorytest == VCF_Source_All ]]; then
	echo "Change directory name and restart"
	exit 1
	fi

for i in *; do
	(if [[ -s $i ]] ; then
        	echo "$i has data" > /dev/null 2>&1
        	else
		echo ""
        	echo "$i is empty.  Fix and restart script"
        	echo ""
		exit 1
	fi
    getbase=`basename "$i"`
    number=`echo $getbase | sed $tbNumberV | sed $tbNumberW`
    echo $number >> list) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
    done
    wait

duplist=`sort list | uniq -d`
rm list
dupNumberSize=`echo $duplist | wc | awk '{print $3}'`
if [ $dupNumberSize -gt 4 ]
then
    echo "These are duplicated VCFs."
    echo "Must remove duplication, and restart script."
    echo $duplist
    exit 1 # Error status
else
    echo "Good! No duplicate VCFs present"
fi
}

#################################################################################

# Looks for defining positions in VCF files.
# If an AC=1 is found at a defined position it is flagged as a posible mixed infection.
# These defining positions must be SNPs found cluster's main branch

function AConeCallPosition () {

positionList=`awk ' { print $2 }' "${DefiningSNPs}" | awk ' NF > 0 '`

echo "AConeCallPosition is running, started -->  `date`"
#echo "*********************************************************************" >> section2
#echo "Possible Mixed Isolates" > section2
#echo "Defining SNPs that are called as AC=1" >> section2
echo "" >> section2
for i in *.vcf; do
(for pos in $positionList; do awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ "^"x"$" ) print FILENAME, "Pos:", $2, "QUAL:", $6, $8 }' $i; done | grep "AC=1;A" | awk 'BEGIN {FS=";"} {print $1, $2}' >> section2) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

#echo "*********************************************************************" >> section2
}

#################################################################################

# This function prepares the filter files.
# awk needs to see a number in the file, so if the file is blank 2 matching numbers are added.  2 numbers because duplicates are kept therefore allowing a number to be pasting into awk when comparing files.

function filterFilespreparation () {

# For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
#python -u /home/tstuber/workspace/scripts/python_scripts/inputXLS.py | sed 's/ u//g' | tr "," "\t" | sed "s/\'//g" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' > ${FilterDirectory}/filterFile.txt

echo "Waiting for filter file creation to complete"
#filterFileCreations
wait
curdr=`pwd`

cd "${FilterDirectory}"

echo "Preparing Filter Files"
for i in *.txt; do
    (getbase=`basename "$i"`
    number=`echo $getbase | sed 's/\(.*\)\..*/\1/'`
    #echo $number
    cat $i | sort | uniq > "$number".num
    if [ $((chromCount)) -eq 1 ]; then
       echo "100000000" >> "$number.num"
       echo "100000000" >> "$number.num"
    elif [ $((chromCount)) -eq 2 ]; then
       echo "chrom1	100000000" >> "$number.num"
       echo "chrom1	100000000" >> "$number.num"
       echo "chrom2	100000000" >> "$number.num"
       echo "chrom2	100000000" >> "$number.num"
    else
        echo "Greater than 2 chromosomes present."
    fi

        rm $i
        mv "$number.num" "$number.txt") &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

cd "${curdr}"
echo "Finished preparing filter files"

}


#################################################################################

# Change SNPs with low QUAL values to N, based on parameter set above in variable settings

function changeLowCalls () {
echo "Changeing low calls, started --> `date`"
for i in *.vcf; do
(base=`basename $i .vcf`; awk -v x=$lowEnd -v y=$highEnd 'BEGIN {OFS="\t"} { if ($6 >= x && $6 <= y) print $1, $2, $3, $4, "N", $6, $7, $8; else print $0 }' $i > ${base}.txt; rm $i; mv ${base}.txt ${base}.vcf) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

}

#################################################################################

function findpositionstofilter () {

echo "`date` --> Finding positions to filter"
# positions have already been filtered via cutting specific positions.
cp filtered_total_pos total_pos
awk '{print $1}' total_pos > prepositionlist
for n  in `cat prepositionlist`; do
	(front=`echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\1/'`
	back=`echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\2/'`
	#echo "front: $front"
	#echo "back: $back"

	positioncount=`awk -v f=$front -v b=$back ' $1 == f && $2 == b {count++} END {print count}' ./*vcf`
	#echo "position count: $positioncount"
	if [ $positioncount -gt 2 ]; then
		#printf "%s\t%s\n" "$front" "$back"
		echo "$n" >> positionlist
	else
		echo $n >> ${d}-DONOT_filtertheseposition.txt
	fi) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

echo "`date` --> Filtering..."
for p in `cat positionlist`; do
	(front=`echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\1/'`
	back=`echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\2/'`
	#echo "front: $front"
	#echo "back: $back"

	maxqual=`awk -v f=$front -v b=$back 'BEGIN{max=0} $1 == f && $2 == b {if ($6>max) max=$6} END {print max}' ./*vcf | sed 's/\..*//'`

	avequal=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $6}' ./*vcf | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`

	maxmap=`awk -v f=$front -v b=$back ' $1 == f && $2 == b {print $8}' ./*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk 'BEGIN{max=0}{if ($1>max) max=$1} END {print max}' | sed 's/\..*//'`

	avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`

	#change maxmap from 52 to 56 2015-09-18
	if [ $maxqual -lt 1300  ] || [ $avequal -lt 800 ] || [ $maxmap -lt 58  ] || [ $avemap -lt 57 ]; then
		echo "maxqual $maxqual" >> filterpositiondetail
		echo "avequal $avequal" >> filterpositiondetail
		echo "maxmap $maxmap" >> filterpositiondetail
		echo "avemap $avemap" >> filterpositiondetail
		echo "position $p" >> filterpositiondetail
		echo ""  >> filterpositiondetail
		echo "$p" >> ${d}-filtertheseposition.txt
	else
		echo "$p" >> ${d}-DONOT_filtertheseposition.txt
		#echo "maxqual $maxqual"
		#echo "maxmap $maxmap"
		#echo "avemap $avemap"
		#echo "position $p"
		#echo ""
	fi) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 10
rm positionlist
rm prepositionlist

rm total_pos

# Filter VCF files
# cat total_pos ${d}-DONOT_filtertheseposition.txt | sort -k1,1n | uniq -d > filtered_total_pos
# fgrep -f filtered_total_pos total_alt > filtered_total_alt
#mv total_alt filtered_total_alt
#rm ${d}-DONOT_filtertheseposition.txt

}

#################################################################################

#   Function: fasta and table creation
function fasta_table () {

# Loop through the directories
directories=`ls`
echo "$directories"
startingdirectory=`pwd`

for d in $directories; do

cd ${startingdirectory}/$d/
dir=`basename $PWD`
echo "Directory:  $dir"

mkdir starting_files
cp *.vcf ./starting_files

	if [ $FilterGroups == yes ]; then
		if [ $((chromCount)) -eq 1 ]; then
			#Mark vcf allowing areas of the genome to be removed from the SNP analysis
			for i in *.vcf; do
				(m=`basename "$i"`; n=`echo $m | sed $dropEXT`
				awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' $i > $i.file
				cat "${FilterDirectory}/$d.txt" $i.file >> $i.catFile
				cat $i.catFile | sort | uniq -d > $i.txt
				pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
				if [ -z $pos ]; then
					#echo "pos is zero... adding a value"
					pos='^1000000000$'
				fi

				awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filtered.vcf

				rm $i.file
				rm $i.catFile
				rm $i.txt
				grep -v "Not_Included" $n.filtered.vcf > $i)  &
				let count+=1
				[[ $((count%NR_CPUS)) -eq 0 ]] && wait
			done
			wait
		else
			#Mark vcf allowing areas of the genome to be removed from the SNP analysis
			for i in *.vcf; do m=`basename "$i"`; n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
				grep '^#' $i > ${i}.header
				grep -v '^#' $i > ${i}.body
				#Mark vcf allowing areas of the genome to be removed from the SNP analysis
				# Iterate through chrom number range
				COUNTER=0
				for c in `cat $dircalled/chroms`; do
					let COUNTER=COUNTER+1
					#echo The counter is $COUNTER
					#echo "********* In $d --> $n working on chromos $c **********"
					awk -v c=$c 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ && $1 == c {print $2}' ${i}.body > $i.filepositions
					awk -v c=$c ' $1 == c {print $2}' ${FilterDirectory}/FilterToAll.txt > $i.positionstofilter
					cat $i.positionstofilter $i.filepositions | sort -k1,1 | uniq -d > $i.foundpositions
					pos=`cat $i.foundpositions | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`

					if [[ -n $pos ]]; then
						echo "pos: $pos" > /dev/null 2>&1
					else
						#echo "string is zero; no findings for pos; giving pos=1"
						pos="^1$"
					fi

					awk -v var1=$c -v var2=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' ${i}.body | grep "$c" > $n.filterchrom${COUNTER}.vcf
				done
				cat ${i}.header $n.filterchrom*.vcf > $n.filtered.vcf
				grep -v "Not_Included" $n.filtered.vcf > $i

				rm ${i}.header
				rm ${i}.body
				rm $i.filepositions
				rm $i.positionstofilter
				rm $n.filterchrom*.vcf
				rm $i.foundpositions

			done
		fi
	mkdir marked_files
	mv *.filtered.vcf ./marked_files
	wait
	sleep 2

	fi

# Make concatemer with the position and REF call.
# Factor in possible multiple chromosomes
# Get rid of duplicates in concatemer and list all the positions and REF calls
for i in *.vcf; do
	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' $i >> concatemer
done

# Get rid of duplicates in concatemer and list all the positions and REF calls
sort -k1,1 < concatemer | uniq > filtered_total_alt
awk '{print $1}' filtered_total_alt > filtered_total_pos

# Count the number of SNPs
totalSNPs=`wc -l  filtered_total_pos`

######################## FILTER FILE CREATOR ###########################
if [ "$cflag" ]; then
	findpositionstofilter
fi
#########################################################################

# Find AC1 positions also found in total_pos
awk '{print $1}' filtered_total_pos > total.list

    for i in *.vcf; do
	(m=`basename "$i"`; n=`echo $m | sed 's/\..*//'`
	# search for AC1 positions
	awk ' $0 !~ /^#/ && $8 ~ /^AC=1/ && $6 > 0 {print $1 "-" $2}' $i > ${n}.list
	# AC1 positions that are being found in this group
	positionsfound=`cat ${n}.list total.list | sort -n | uniq -d`
	countfind=`echo $positionsfound | wc -w`
	#echo "positonsfound: $positionsfound  countfind: $countfind"
	rm ${n}.list
	if [[ -z $positionsfound ]]; then
		positionsfound="No positions found"
	fi
	
	if [[ $countfind -gt 2  ]]; then
	searchname=`echo $n | sed 's/_.*//'`

        	if [[  $argUsed == para ]]; then
            	unmappedContigs=`grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/data/${searchname}/bwamem-gatk/qualityvalues/*stats.txt`
        	elif [[  $argUsed == bovis ]]; then
            	unmappedContigs=`grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script1/${searchname}/bwamem-gatk/qualityvalues/*stats.txt`
        	else
            contigMessage="possibly set a new contig path at script line: $LINENO"
        	fi

        	if [[ -z $unmappedContigs ]]; then
			unmappedContigs="Contig counts not available"
		fi

	echo "$d" >> $d-AC1postions.txt
	echo "" >> $d-AC1postions.txt

	echo "$d Sample: $n  AC1 findings: $countfind  $unmappedContigs $contigMessage" > delete
	tr -d "\n" < delete >> $d-AC1postions.txt
	echo "" >> $d-AC1postions.txt
	tr -d "\n" < delete >> ${fulDir}/emailAC1counts.txt
        echo "" >> ${fulDir}/emailAC1counts.txt
	
	for p in `echo $positionsfound`; do
            position=`echo $p | sed 's/chrom[0-9]*-//'`
            awk -v p="$position" '$2 == p {print $0}' $i >> $d-AC1postions.txt
        done

    fi)  &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
    done
    wait
rm total.list
rm delete
# Count the number of SNPs

totalSNPs=`grep -c ".*" filtered_total_pos`
echo "$d total SNPs: $totalSNPs" >> ../../section4

echo "***Creating normalized vcf using AC2, QUAL > $QUAL"
# Grab the name of the vcf file

#########################################################################
# Count the number of SNPs
filteredSNPs=`wc -l filtered_total_pos`
echo "Total SNPs after filtering $filteredSNPs"


for i in *.vcf; do
	(n=${i%.vcf}
	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' $i > $n.allsnps_alt
	#get SNPs of interest
	fgrep -f filtered_total_pos $n.allsnps_alt > $n.targetsnps_alt
	#if SNP not found in sample default call to reference, normalize.
	cat $n.targetsnps_alt filtered_total_alt | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1,1 > $n.filteredsnps_alt

	# If position has zero map quality change alt call to -
	# get positions being used
	awk '{print $1}' $n.filteredsnps_alt > $n.filteredsnps_pos
	# Get zero coverage positions.
	awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' ${i} > ${n}.zeropositions

	# if duplicate then zero mapped position found for sample
	cat $n.filteredsnps_pos ${n}.zeropositions | sort | uniq -d | awk '{print $1, "-"}' > ${n}.zerotomerge_alt #the - makes it and alt file

	#if zero positions found merge them to the SNPs found
	if [ -s ${n}.zerotomerge_alt ]; then
		# merge zero updates to SNP file
		cat ${n}.zerotomerge_alt $n.filteredsnps_alt | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > ${n}.zerofilteredsnps_alt
		#echo "***Found zero postions: $n"
		rm $n.filteredsnps_alt
	else
		#echo "no zero positions found for $n"
		mv $n.filteredsnps_alt ${n}.zerofilteredsnps_alt
	fi

	rm $n.allsnps_alt
	rm $n.filteredsnps_pos
	rm $n.targetsnps_alt
	rm ${n}.zeropositions
	rm ${n}.zerotomerge_alt)  &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 5
wait

echo "`date` --> Finding parsimony informative positions"
# Capture only positions that have more than one SNP type called at a position
cat *zerofilteredsnps_alt | sort -nk1,1 | uniq | awk '{print $1}' | uniq -d > parsimony_informative
# This removes calls that are the same for all isolates being analyzed

# If many SNPs fgrep may not do much and be slow
fgrep -f parsimony_informative filtered_total_alt | sort -k1,1n > parsimony_filtered_total_alt
awk '{print $1}' parsimony_filtered_total_alt > parsimony_filtered_total_pos

# Create table and fasta
awk '{print $1}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt
awk '{print $2}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt

for i in *zerofilteredsnps_alt; do
	(m=`basename "$i"`; n=`echo $m | sed 's/\..*//'`

	fgrep -f parsimony_filtered_total_pos $i | sort -k1,1n > $n.pretod

	##############################################################
	# Change AC1s to IUPAC

	# get positions being used
	awk '{print $1}' ${n}.pretod > ${n}.usedpostions
	# get AC1 positions and iupac calls  that were changed to iupac
	awk -v Q="$QUAL" ' $0 !~ /#/ && $6 > Q && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' ${i%zerofilteredsnps_alt}vcf > ${n}.ac
	# get just positions of those AC1 grabbed above
	awk '{print $1}' ${n}.ac > ${n}.acpositions
	# AC duplicate positions will need to be kept
	cat ${n}.usedpostions ${n}.acpositions | sort | uniq -d > ${n}.actokeep
	# get AC1 position with iupac, these are only positions already in the pretod

	if [ -s ${n}.actokeep ]; then
		fgrep -f ${n}.actokeep ${n}.ac > ${n}.actomerge
		# merge iupac updates to filledcut
		cat ${n}.actomerge $n.pretod | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > $n.tod
		rm ${n}.pretod
		rm ${n}.actomerge
	else
		#echo "else done"
		mv $n.pretod $n.tod
	fi
	rm ${n}.usedpostions
	rm ${n}.ac
	rm ${n}.acpositions
	rm ${n}.actokeep
	##############################################################

	awk '{print $2}' $n.tod | tr -d [:space:] | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > $n.fas
	# Add each isolate to the table
	awk '{print $2}' $n.tod | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt ) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

wait
sleep 5

#Create root sequence
awk '{print $2}' parsimony_filtered_total_alt > root
cat root | tr -cd "[:print:]" | sed "s/^/>root;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > root.fas
echo "" >> root.fas

totalSNPs=`grep -c ".*" parsimony_filtered_total_pos`
echo "Total informative SNPs: $totalSNPs"

#Clean-up
rm concatemer
rm *.tod
mkdir fasta
mv *.fas ./fasta
rm root
rm *vcf
rm filtered_total_alt
rm filtered_total_pos
rm parsimony_filtered_total_alt
rm parsimony_filtered_total_pos
rm parsimony_informative
rm *zerofilteredsnps_alt

cp /home/shared/Table_Template.xlsx ./${d}-Table_Template.xlsx
done
}
#****************************************************************
function alignTable () {

# Beginning in fasta folder
echo "`date` --> RAxML started $d"

awk '{print $0}' *.fas | sed '/root/{N;d;}' >> fastaGroup.txt
awk '{print $0}' *.fas >> RAxMLfastaGroup.txt

raxmlHPC-SSE3 -s RAxMLfastaGroup.txt -n ${d} -m GTRCAT -p 12345 &>/dev/null && nw_reroot RAxML_bestTree.${d} root | nw_display -s -w 1000 -v 20 -b 'opacity:0' -i 'font-size:8' -l 'font-family:serif;font-style:italic' -d 'stroke-width:2;stroke:blue' - > ../${d}-tree.svg; nw_reroot RAxML_bestTree.${d} root > tableinput.${d}; nw_reroot RAxML_bestTree.${d} root > rooted_RAxML_bestTree.${d}; mv rooted_RAxML_bestTree.${d} RAxML_bestTree.${d}
wait
rm RAxML_parsimonyTree*
for i in RAxML*Tree*; do mv $i ../${i}.tre; done

tr ":" "\n" < tableinput.${d} | tr "," "\n" | sed 's/(//g' | sed 's/)//g' | grep -v "\.[0-9]*" | grep -v "root" > cleanedAlignment.txt

awk 'NR==FNR{o[FNR]=$1; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' cleanedAlignment.txt ../$d.table.txt > joined.txt
grep "reference" ../$d.table.txt > references.txt
cat references.txt joined.txt >> joined2.txt
mv joined2.txt ../$d.sortedTable.txt

rm joined.txt
rm references.txt

cd ..

#get the number of columns
columnCount=`awk '$0 !~ /^$/ {print $0}' *sortedTable.txt | awk '{print NF-1; exit}'`

#number=`jot - 1 $columnCount`
number=`seq $columnCount`
#echo "Numbers in list: $number"

#get the reference row
cat *sortedTable.txt | sed -n 2p | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > referenceRow
cat referenceRow > out2.txt
cat *sortedTable.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > table
#remove first column from *sortedTable.txt

echo "countDif" > countOutput.txt
echo "countFirst" > firstOutput.txt

#iterate numbers up to number of columns
for n in $number; do
	#use number to grab character in reference row i.e. C
	letter=`awk -v x=$n '{print $x}' referenceRow`
	#echo "Number: $n"
	#echo "Letter: $letter"
	awk -v var1="$n" '{print $var1}' table > column
	grep "$letter" column | wc -l >> countOutput.txt
	sed '1,2d' column > column2
	cat column2 | awk -v var2=${letter} ' $0 !~ var2 {print NR; exit}' >> firstOutput.txt
done

#/var2/ {print NR; exit}' column
#Clear table2.txt
echo "" > table2.txt
#Add a \n to the end of file

cat *sortedTable.txt >> table2.txt
#sed '$a\' *sortedTable.txt >> table2.txt

#Prepare count line
cat countOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyline.txt
#Create readytable.txt with count line.
cat table2.txt readyline.txt > orginizedTable2.txt
grep -v "^$" orginizedTable2.txt > orginizedTable3.txt

#Prepare firstOut line
cat firstOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyFirstOut.txt
cat orginizedTable3.txt readyFirstOut.txt > orginizedTable4.txt

rm referenceRow
rm out2.txt
rm table
rm countOutput.txt
rm table2.txt
rm readyline.txt
rm orginizedTable2.txt

awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable4.txt > orginizedTable5.txt
awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable5.txt > orginizedTable6.txt

#Transpose
awk '{
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str" "a[i,j];
}
print str
}
}' orginizedTable6.txt > orginizedTable7.txt

#Orgainize file based on 1st 2 columns
sort -n -k1 orginizedTable7.txt | sort -n -k2 > orginizedTable8.txt

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable8.txt > orginizedTable9.txt
awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' orginizedTable9.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > orginizedTable10.txt

#Transpose back
awk '{
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str" "a[i,j];
}
print str
}
}' orginizedTable10.txt > orginizedTable11.txt

c=`basename $PWD`

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable11.txt > $c.organizedTable.txt

rm orginizedTable3.txt
rm orginizedTable4.txt
rm orginizedTable5.txt
rm orginizedTable6.txt
rm orginizedTable7.txt
rm orginizedTable8.txt
rm orginizedTable9.txt
rm orginizedTable10.txt
rm orginizedTable11.txt
rm column
rm column2
rm readyFirstOut.txt
rm firstOutput.txt
echo "Adding map qualities..."

# Add map qualities to sorted table

# Get just the position.  The chromosome must be removed
awk ' NR == 1 {print $0}' $d.sortedTable.txt | tr "\t" "\n" | sed "1d" | awk '{print NR, $0}' > $d-positions

echo "map-quality map-quality" > quality.txt
echo "`date` --> Sorted table map quality gathering for $c"
while read p; do
	(rownumber=`echo $p | awk '{print $1}'`
	front=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/'`
	back=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/'`
	#echo "rownumber: $rownumber"
	#echo "front: $front"
	#echo "back: $back"
	avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./starting_files/*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
	echo "$rownumber $avemap" >> quality.txt) &
	CPU_NR=$(mpstat | grep -A 5 "%idle" | tail -n 1 | awk -F " " '{print  64 * (0.01 * $12) - 15}'a | sed 's/\..*//')
	let count+=1
	[[ $((count%CPU_NR)) -eq 0 ]] && wait
	done < $d-positions

sort -nk1,1 < quality.txt | awk '{print $2}' | tr "\n" "\t" > qualitytransposed.txt

cat $d.sortedTable.txt qualitytransposed.txt | grep -v '^$' > $d-mapquality-sortedtable.txt
mv $d-mapquality-sortedtable.txt $d.sortedTable.txt

# Add map qualities to organized table

awk ' NR == 1 {print $0}' $c.organizedTable.txt | tr "\t" "\n" | sed "1d" | awk '{print NR, $0}' > $d-positions

echo "map-quality map-quality" > quality.txt
echo "`date` --> Organized table map quality gathering for $c"
while read p; do
	(rownumber=`echo $p | awk '{print $1}'`
	front=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/'`
	back=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/'`
	#echo "rownumber: $rownumber"
	#echo "front: $front"
	#echo "back: $back"
	avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./starting_files/*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
	echo "$rownumber $avemap" >> quality.txt) &
        CPU_NR=$(mpstat | grep -A 5 "%idle" | tail -n 1 | awk -F " " '{print  64 * (0.01 * $12) - 15}'a | sed 's/\..*//')
	let count+=1
	[[ $((count%CPU_NR)) -eq 0 ]] && wait
	done < $d-positions

sort -nk1,1 < quality.txt | awk '{print $2}' | tr "\n" "\t" > qualitytransposed.txt

cat $c.organizedTable.txt qualitytransposed.txt | grep -v '^$' > $d-mapquality-orgainizedtable.txt
mv $d-mapquality-orgainizedtable.txt $c.organizedTable.txt

rm quality.txt
rm qualitytransposed.txt
rm $d-positions

}

#################################################################################
#################################################################################
#################################################################################
###################################### START ####################################
#################################################################################
#################################################################################
#################################################################################

# Clean the tag file that has been exported to Desktop
chmod 777 ${genotypingcodes}  
cat ${genotypingcodes} | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' | sed 's/MN_Wildlife_Deer_//' > preparedTags.txt

#clean_tag.sh $genotypingcodes
####################
# Clean the genotyping codes used for naming output
sed 's/\*//g' < preparedTags.txt | sed 's/(/_/g' | sed 's/)/_/g' | sed 's/ /_/g' | sed 's/-_/_/g' | sed 's/\?//g' | sed 's/_-/_/g' | sed 's/,/_/g' | sed 's#/#_#g' | sed 's#\\#_#g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/-$//g' | sed 's/_$//g' |awk 'BEGIN {OFS="\t"}{gsub("_$","",$1)}1' > outfile
rm preparedTags.txt

cat ${genotypingcodes} | tr '\r' '\n' | grep "Yes" | sed 's/_.*//' >> elite
echo "Only samples in this file will be ran when elite is used as the secound argument" >> elite

####################

# Test for duplicate VCFs
testDuplicates
wait
#Prepare Filter files.
#filterFilespreparation
wait
#Test for match coverage file
#checkMatchingCoverageFile

#copy the original vcfs to /starting_files
mkdir starting_files
mv *.* ./starting_files

# If bovis are ran default will only run with files check "misc" in FileMaker
# Untitled.tab exported from FileMaker must contain "isolate names" followed by "Misc".

	if [ "$eflag" ]; then
        echo "Only analyzing elite files"

        for i in `cat elite`; do
        name=`ls starting_files | grep $i`
        cp ./starting_files/$name ./
        done

        for i in `find ./starting_files/ -mtime -10`; do
        cp $i ./
        done

    else
        echo "all samples will be ran"
        cp ./starting_files/* ./
	fi

rm elite

#Remove possible "## in vcf headers
echo 'Removing possible "## in vcf headers'
sed -i 's/^"##/##/' *vcf

#################################################################################

# Count the number of chromosomes used in the reference when VCFs were made.
#singleFile=`ls *.vcf | head -1`
echo "Counting the number of chromosomes in first 100 samples, started -->  `date`"
chromCount=`awk ' $0 !~ /^#/ {print $1}' $(ls *vcf | head -100) | sort | uniq -d | awk 'END {print NR}'`
echo "The number of chromosomes/segments seen in VCF: $chromCount"
awk ' $0 !~ /^#/ {print $1}' $(ls *vcf | head -100) | sort | uniq -d > chroms
echo "These are the chromosomes/segments found:"
cat chroms

#################################################################################

# Remove selected isolates from comparison
# This is optional, and should be turned on or off based on laboratories preference
removeIsolates

############################### Rename files ###############################

for i in *.txt; do
    mv $i ${i%.txt}.vcf
done
for i in *.vcf; do
    echo "******************** Naming convention ********************"
    echo "Original File: $i"
    base=`basename "$i"`
    searchName=`echo $base | sed $tbNumberV | sed $tbNumberW | sed 's/V//'`
    echo "searchName: $searchName"
    # Direct script to text file containing a list of the correct labels to use.
    # The file must be a txt file.
    p=`grep "$searchName" "outfile"`
    echo "This is what was found in tag file: $p"
    newName=`echo $p | awk '{print $1}' | tr -d "[:space:]"` # Captured the new name
    n=`echo $base | sed $tbNumberV | sed $tbNumberW`
    noExtention=`echo $base | sed $dropEXT`
    VALtest=`echo $i | grep "VAL"`
#    echo "VALtest: $VALtest"
#h=`echo ${i%-AZ}`; g=`echo ${h%-Broad}`; echo $g
    #Check if a name was found in the tag file.  If no name was found, keep original name, make note in log and cp file to unnamed folder.
    if [[ -z "$p" ]]; then # new name was NOT found
        if [[ -z "$VALtest" ]]; then
            name=$searchName
            echo "n is $n"
            echo "$name" >> section1
            mkdir -p FilesNotRenamed
            cp $i ./FilesNotRenamed
            mv $i ${name}.vcf
#            echo "A"
        else
            name=${searchName}-Val
            mv $i ${name}.vcf
#            echo "B"
        fi
    else # New name WAS found
        if [[ -z "$VALtest" ]]; then
            name=$newName
            mv $i ${name}.vcf
#            echo "C"
        else
            name=${newName}-Val
            echo "newName is $name"
            mv $i ${name}.vcf
#            echo "D"
        fi
    fi
done
rm outfile
##################### Start: Make Files Unix Compatiable #####################

#Fix validated (VAL) vcf files.  This is used in vcftofasta scripts to prepare validated vcf files opened and saved in Excel.
#Create list of isolates containing "VAL"
#Do NOT make this a child process.  It messes changing column 1 to chrom

echo "#############################"
echo "Making Files Unix Compatiable"
for v in *.vcf; do
    (dos2unix $v > /dev/null 2>&1 #Fixes files opened and saved in Excel
    cat $v | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' > $v.temp
    mv $v.temp $v) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

########################################################################

printf "%s\t%s\t%s\t%s\n" "TB Number" "Group" "Subgroup" "Clade" > FileMakerGroupImport.txt

AConeCallPosition
wait

# Change low QUAL SNPs to N, see set variables
changeLowCalls
wait

######################## Change AC1s to IUPAC ########################

echo "Changing AC=1 to IUPAC, started -->  `date`"

for i in *vcf; do

(awk '
BEGIN { OFS = "\t"}

{ if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AG/ )
	 	print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CT/ )
	 	print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GC/ )
	 	print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AT/ )
	 	print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GT/ )
	 	print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AC/ )
	 	print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GA/ )
	 	print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TC/ )
	 	print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CG/ )
	 	print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TA/ )
	 	print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TG/ )
	 	print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10	 	
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CA/ )
	 	print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10	 	
else
	print $0	 
}' $i > ${i%vcf}temp
mv ${i%vcf}temp $i) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
######################## Mark Files and Remove Marked Regions ########################

if [ $FilterAllVCFs == yes ]; then
echo "`date` --> Marking all VCFs and removing filtering region"
	# Label filter field for positions to be filtered in all VCFs
        if [ $((chromCount)) -eq 1 ]; then
        for i in *.vcf; do
        (m=`basename "$i"`; n=`echo $m | sed $dropEXT`
        # Get usable positions in the VCF
        awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' $i > $i.file
        # Combine with positions that will be filtered
        cat "${FilterDirectory}/FilterToAll.txt" $i.file >> $i.catFile
        
	# Output any duplicate positions, aka decreasing positions to be marked and used by awk
        cat $i.catFile | sort | uniq -d > $i.txt
        # preparing postions
        pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`

	# Making a vcf with positions marked that should not be included based on filter file
        awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filtered.vcf
        # Clean up

        rm $i.file; rm $i.catFile; rm $i.txt
	# Removed positions and clobber origingal vcf
	grep -v "Not_Included" $n.filtered.vcf > $i) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
	done
	wait

	elif [ $((chromCount)) -gt 1 ]; then
		#echo "multiple chromosomes"
		for i in *.vcf; do m=`basename "$i"`; n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
		grep '^#' $i > ${i}.header
		grep -v '^#' $i > ${i}.body
		#Mark vcf allowing areas of the genome to be removed from the SNP analysis
		# Iterate through chrom number range
		COUNTER=0
		for c in `cat chroms`; do
			let COUNTER=COUNTER+1
			#echo The counter is $COUNTER
			#echo "********* In all_vcfs $n working on chromos $c **********"
			awk -v c=$c 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ && $1 == c {print $2}' ${i}.body > $i.filepositions
			awk -v c=$c ' $1 == c {print $2}' ${FilterDirectory}/FilterToAll.txt > $i.positionstofilter
			cat $i.positionstofilter $i.filepositions | sort -k1,1 | uniq -d > $i.foundpositions
			pos=`cat $i.foundpositions | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
			if [[ -n $pos ]]; then
				echo "pos: $pos" > /dev/null 2>&1
			else
				#echo "string is zero; no findings for pos; giving pos=1"
				pos="^1$"
				#echo $pos
			fi
			awk -v var1=$c -v var2=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' ${i}.body | grep "$c" > $n.filterchrom${COUNTER}.vcf
		done
			cat ${i}.header $n.filterchrom*.vcf > $n.filtered.vcf
			grep -v "Not_Included" $n.filtered.vcf > $i
			rm ${i}.header
			rm ${i}.body
			rm $i.filepositions
			rm $i.positionstofilter
			rm $n.filterchrom*.vcf
			rm $i.foundpositions
		done
		else
		echo "Check chromosome count numbers at line $LINENO.  Exiting script."
		exit 1
	fi

wait
sleep 2

mkdir marked_files
mv *.filtered.vcf ./marked_files

    else
    echo "***All VCF filtering was NOT done."
    echo "***All VCF filtering was NOT done." >> section5
fi

#################### Categorize VCFs into Groups, Subgroups and Clades #####################

# Print header
#    printf "%10s %10s %10s %10s\n" "Isolate" "Group" "Subgroup" "Clade" >> log

echo "" > section3
echo "NAME GROUP SUBGROUP CLADE" >> section3
echo "" >> section3

for i in *.vcf; do

	# If there is one chromosome present just get the position.  If multiple chromosomes are present than the chromsome identification needs to be identified.  The filter file needs to sync with this chromosome identification.  If multiple chromosomes the filter file will be kept as a text file.  If a single chromosome an linked Excel file can be used.
	if [ $((chromCount)) -eq 1 ]; then
		# Get quality positions in VCF
		awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/{print $2}' $i > quality-${i%.vcf}
		else
		# Get quality positions in VCF and include chromosome identification
		awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2}' $i > quality-${i%.vcf}
	fi

	echo "quality-${i%.vcf}:"

	##----Group

	# If a group number matches a quality position in the VCF (formatedpos) then print the position
	grep "Group" "${DefiningSNPs}" > groupsnps

	awk 'NR==FNR{a[$0];next}$2 in a' quality-${i%.vcf} groupsnps | awk '{print $1}' > group-foundpositions-${i%.vcf}

	echo "This is the Group Numbers: `cat group-foundpositions-${i%.vcf}`"

	# Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
	sizeGroup=`wc -l group-foundpositions-${i%.vcf} | awk '{print $1}'`

	# Loop through the number of groups positions found
	loops=`cat group-foundpositions-${i%.vcf}`

	if [ $sizeGroup -lt 1 ]; then # There was not a position found that places VCF into group
		echo "$i Grp not found" >> section3
		echo "$i was not assigned a Group"
		elif [ $sizeGroup -gt 1 ]; then
		echo "$i has multiple groups" >> section3
		echo "$i has multiple groups"
		for l in $loops; do
			echo "making group $i"
			mkdir -p Group-$l #Make groupNumber folder if one does not exist.
			cp $i ./Group-$l/ #Then copy to each folder
		done
		else
		echo "Just one group"
		mkdir -p Group-$loops #Make groupNumber folder if one does not exist.
		cp $i ./Group-$loops/ #Then copy to each folder

	fi

	##----Subgroup

	# If a group number matches a quality position in the VCF (formatedpos) then print the position
	grep "Subgroup" "${DefiningSNPs}" > subgroupsnps

	awk 'NR==FNR{a[$0];next}$2 in a' quality-${i%.vcf} subgroupsnps | awk '{print $1}' > subgroup-foundpositions-${i%.vcf}

	echo "This is the Subgroup Numbers: `cat subgroup-foundpositions-${i%.vcf}`"

	# Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
	sizeGroup=`wc -l subgroup-foundpositions-${i%.vcf} | awk '{print $1}'`

	# Loop through the number of groups positions found
	loops=`cat subgroup-foundpositions-${i%.vcf}`

	if [ $sizeGroup -lt 1 ]; then # There was not a position found that places VCF into group
		echo "$i was not assigned a Subgroup"
		elif [ $sizeGroup -gt 1 ]; then
		echo "$i has multiple subgroups" >> section3
		echo "$i has multiple subgroups"
		for l in $loops; do
			echo "making subgroup $i"
			mkdir -p Subgroup-$l #Make groupNumber folder if one does not exist.
			cp $i ./Subgroup-$l/ #Then copy to each folder
		done
		else
		echo "Just one Subgroup"
		mkdir -p Subgroup-$loops #Make groupNumber folder if one does not exist.
		cp $i ./Subgroup-$loops/ #Then copy to each folder

	fi

	##----Clade

	# If a group number matches a quality position in the VCF (formatedpos) then print the position
	grep "Clade" "${DefiningSNPs}" > cladesnps

	awk 'NR==FNR{a[$0];next}$2 in a' quality-${i%.vcf} cladesnps | awk '{print $1}' > clade-foundpositions-${i%.vcf}

	echo "This is the Clade Numbers: `cat clade-foundpositions-${i%.vcf}`"

	# Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
	sizeGroup=`wc -l clade-foundpositions-${i%.vcf} | awk '{print $1}'`

	# Loop through the number of groups positions found
	loops=`cat clade-foundpositions-${i%.vcf}`

	if [ $sizeGroup -lt 1 ]; then # There was not a position found that places VCF into group
		echo "$i was not assigned a Clade"
		elif [ $sizeGroup -gt 1 ]; then
		echo "$i has multiple clades" >> section3
		echo "$i has multiple clades"
		for l in $loops; do
			echo "making clade $i"
			mkdir -p Clade-$l #Make groupNumber folder if one does not exist.
			cp $i ./Clade-$l/ #Then copy to each folder
		done
		else
		echo "Just one clade"
		mkdir -p Clade-$loops #Make groupNumber folder if one does not exist.
		cp $i ./Clade-$loops/ #Then copy to each folder
	fi
	echo "${i%.vcf} $(cat group-foundpositions-${i%.vcf} subgroup-foundpositions-${i%.vcf} clade-foundpositions-${i%.vcf})" | tr "\n" "\t" >> section3
	echo "" >> section3

	echo ""
	rm quality-${i%.vcf}
	rm groupsnps
	rm subgroupsnps
	rm cladesnps
	rm *foundpositions-${i%.vcf}
	######

	mkdir -p all_vcfs #Make all_vcfs folder if one does not exist.
	mv $i ./all_vcfs/

done

################### Organize folders #####################

mkdir all_groups
mv ./Group-*/ ./all_groups
mkdir all_subgroups
mv ./Subgroup*/ ./all_subgroups/
mkdir all_clades
mv ./Clade*/ ./all_clades/

##################### Start: All vcf folder #####################
cd ./all_vcfs/
d="all_vcfs"

# Make concatemer with the position and REF call.
# Factor in possible multiple chromosomes
# Get rid of duplicates in concatemer and list all the positions and REF calls
echo "`date` --> Gathering SNP positions"

for i in *.vcf; do
	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' $i >> concatemer
done

# Get rid of duplicates in concatemer and list all the positions and REF calls
sort -k1,1 < concatemer | uniq > filtered_total_alt
awk '{print $1}' filtered_total_alt > filtered_total_pos

# Count the number of SNPs
totalSNPs=`wc -l  filtered_total_pos`
echo "Total filtered SNPs: $totalSNPs"

for i in *.vcf; do
	(n=${i%.vcf}
	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' $i > $n.allsnps_alt
	#get SNPs of interest
	fgrep -f filtered_total_pos $n.allsnps_alt > $n.targetsnps_alt
	#if SNP not found in sample default call to reference, normalize.
	cat $n.targetsnps_alt filtered_total_alt | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1,1 > $n.filteredsnps_alt

	# If position has zero map quality change alt call to -
	# get positions being used
	awk '{print $1}' $n.filteredsnps_alt > $n.filteredsnps_pos
	# Get zero coverage positions.
	awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' ${i} > ${n}.zeropositions

	# if duplicate then zero mapped position found for sample
	cat $n.filteredsnps_pos ${n}.zeropositions | sort | uniq -d | awk '{print $1, "-"}' > ${n}.zerotomerge_alt #the - makes it and alt file

	#if zero positions found merge them to the SNPs found
	if [ -s ${n}.zerotomerge_alt ]; then
		# merge zero updates to SNP file
		cat ${n}.zerotomerge_alt $n.filteredsnps_alt | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > ${n}.zerofilteredsnps_alt
		#echo "***Found zero postions: $n"
		rm $n.filteredsnps_alt
	else
		#echo "no zero positions found for $n"
		mv $n.filteredsnps_alt ${n}.zerofilteredsnps_alt
	fi

	rm $n.allsnps_alt
	rm $n.filteredsnps_pos
	rm $n.targetsnps_alt
	rm ${n}.zeropositions
	rm ${n}.zerotomerge_alt)  &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 5
wait

echo "`date` --> Finding parsimony informative positions"
# Capture only positions that have more than one SNP type called at a position
cat *zerofilteredsnps_alt | sort -nk1,1 | uniq | awk '{print $1}' | uniq -d > parsimony_informative
# This removes calls that are the same for all isolates being analyzed
# If many SNPs fgrep may not do much and be slow
fgrep -f parsimony_informative filtered_total_alt | sort -k1,1n > parsimony_filtered_total_alt
awk '{print $1}' parsimony_filtered_total_alt > parsimony_filtered_total_pos

######################## FILTER FILE CREATOR ###########################
if [ "$cflag" ]; then
	d="all_vcfs"
	findpositionstofilter
fi
#########################################################################

# Create table and fasta
awk '{print $1}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt
awk '{print $2}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt

for i in *zerofilteredsnps_alt; do
	(m=`basename "$i"`; n=`echo $m | sed 's/\..*//'`

	fgrep -f parsimony_filtered_total_pos $i | sort -k1,1n > $n.pretod

	##############################################################
	# Change AC1s to IUPAC

	# get positions being used
	awk '{print $1}' ${n}.pretod > ${n}.usedpostions
	# get AC1 positions and iupac calls  that were changed to iupac
	awk -v Q="$QUAL" ' $0 !~ /#/ && $6 > Q && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' ${i%zerofilteredsnps_alt}vcf > ${n}.ac
	# get just positions of those AC1 grabbed above
	awk '{print $1}' ${n}.ac > ${n}.acpositions
	# AC duplicate positions will need to be kept
	cat ${n}.usedpostions ${n}.acpositions | sort | uniq -d > ${n}.actokeep
	# get AC1 position with iupac, these are only positions already in the pretod

	if [ -s ${n}.actokeep ]; then
		fgrep -f ${n}.actokeep ${n}.ac > ${n}.actomerge
		# merge iupac updates to filledcut
		cat ${n}.actomerge $n.pretod | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > $n.tod
		rm ${n}.pretod
		rm ${n}.actomerge
	else
		#echo "else done"
		mv $n.pretod $n.tod
	fi
	rm ${n}.usedpostions
	rm ${n}.ac
	rm ${n}.acpositions
	rm ${n}.actokeep
	##############################################################

	awk '{print $2}' $n.tod | tr -d [:space:] | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > $n.fas
	# Add each isolate to the table
	awk '{print $2}' $n.tod | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt ) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

wait
sleep 5

#Create root sequence
awk '{print $2}' parsimony_filtered_total_alt > root
cat root | tr -cd "[:print:]" | sed "s/^/>root;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > root.fas
echo "" >> root.fas

totalSNPs=`grep -c ".*" parsimony_filtered_total_pos`
echo "Total informative SNPs: $totalSNPs"

#Clean-up
rm concatemer
rm *.tod
mkdir fasta
mv *.fas ./fasta
rm root
rm *vcf
rm filtered_total_alt
rm filtered_total_pos
rm parsimony_filtered_total_alt
rm parsimony_filtered_total_pos
rm parsimony_informative
rm *zerofilteredsnps_alt

if [ "$eflag" -o "$aflag" ]; then
	d="all_vcfs"
        cd ./fasta
        alignTable
else
	echo "Tree not ran for all_vcfs"
fi

##################### End: All vcf folder #####################

#echo "***************************************************"
#echo "***************** STARTING Groups *****************"
#echo "***************************************************"
# Change directory to all_groups
cd ${fulDir}/all_groups
fasta_table &  

#echo "***************************************************"
#echo "**************** STARTING SUBGROUPS ***************"
#echo "***************************************************"
# Change directory to all_subgroups
cd ${fulDir}/all_subgroups
fasta_table &

#echo "***************************************************"
#echo "***************** STARTING CLADES *****************"
#echo "***************************************************"
# Change directory to all_clades
cd ${fulDir}/all_clades
fasta_table &
wait
echo "At line $LINENO, sleeping 5 second"; sleep 5s
cd ${fulDir}
$PWD
cp ${DefiningSNPs} ./
cp /home/shared/Table_Template.xlsx ./
cp "$0" "$PWD"

#echo "***************************************************"
#echo "********** STARTING all_clades Alignment **********"
#echo "***************************************************"
cd ${fulDir}/all_clades

workingdir=`basename $PWD`

if [ $workingdir == all_clades ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_clades/${d}/fasta
    #echo "****************************************************"
    #echo "************* Orginizing Table: $d *****************"
    #echo "****************************************************"
	alignTable & 

    pwd
done

else
echo "*** $workingdir not found ***"
fi

#echo "***************************************************"
#echo "********** STARTING all_groups Alignment **********"
#echo "***************************************************"
cd ${fulDir}/all_groups

workingdir=`basename $PWD`

if [ $workingdir == all_groups ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_groups/${d}/fasta
    #echo "****************************************************"
    #echo "************* Orginizing Table: $d *****************"
    #echo "****************************************************"
    alignTable &
pwd
done
else
echo "*** $workingdir not found ***"
fi

#echo "***************************************************"
#echo "******** STARTING all_subgroups Alignment *********"
#echo "***************************************************"
cd ${fulDir}/all_subgroups

workingdir=`basename $PWD`

if [ $workingdir == all_subgroups ]
then

directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_subgroups/$d/fasta
    #echo "****************************************************"
    #echo "************* Orginizing Table: $d *****************"
    #echo "****************************************************"
    alignTable &
pwd
done
else
echo "*** $workingdir not found ***"
fi
wait
pwd
cd ${fulDir}

column section1 > csection1
sort -nr < section4 > ssection4

echo "End Time:  `date`" >> sectiontime
endtime=`date +%s`
runtime=$((endtime-starttime))
#totaltime=`date -u -d @${runtime} +"%T"`
echo "Run time: $runtime seconds" >> sectiontime

cat sectiontime >  log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "" >> log.txt
cat section5 >> log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "" >> log.txt
echo "These files did not get renamed:" >> log.txt
cat csection1 >> log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "" >> log.txt
echo "Possible Mixed Isolates" >> log.txt
echo "Defining SNPs called AC=1" >> log.txt
cat section2 >> log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "" >> log.txt
cat section3 >> log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "SNP counts::" >> log.txt
cat ssection4 >> log.txt
echo "" >> log.txt
echo "****************************************************" >> log.txt
echo "AC1 called SNPs"
cat ${fulDir}/emailAC1counts.txt | sort -nk1,1 >> log.txt

echo "<html>" > email_log.html
echo "<Body>" >> email_log.html
awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' sectiontime >  email_log.html
echo "****************************************************" >> email_log.html
echo "" >> email_log.html
awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' section5 >> email_log.html
echo "" >> email_log.html
echo "****************************************************" >> email_log.html
echo "" >> email_log.html
echo "<p> These files did not get renamed: </p>" >> email_log.html
awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' csection1 >> email_log.html
echo "" >> email_log.html
echo "****************************************************" >> email_log.html
echo "" >> email_log.html
echo "<p> Possible Mixed Isolates, Defining SNPs called AC=1 </p>" >> email_log.html
awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' section2 >> email_log.html
echo "" >> email_log.html
echo "****************************************************" >> email_log.html
echo "" >> email_log.html
awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' section3 >> email_log.html
echo "" >> email_log.html
echo "****************************************************" >> email_log.html
echo "" >> email_log.html
echo "<p> SNP counts: </p>" >> email_log.html
awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' ssection4 >> email_log.html
echo "" >> email_log.html
echo "****************************************************" >> email_log.html
echo "<p> AC1 called SNPs: </p>" >> email_log.html
awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' ${fulDir}/emailAC1counts.txt >> email_log.html
echo "</Body>" >> email_log.html
echo "</html>" >> email_log.html

rm section1
rm section2
rm section3
rm section4
rm section5
rm sectiontime
rm ssection4
rm csection1
#rm -r ${FilterDirectory}

echo "Copying to ${bioinfoVCF}"
cp -r $PWD ${bioinfoVCF}
fileName=`basename $0`

if [ "$mflag" ]; then
    email_list="Tod.P.Stuber@aphis.usda.gov"
	echo "vcftofasta.sh completed" > mytempfile; cat mytempfile | mutt -s "vcftofasta.sh completed subject" -a email_log.html -- $email_list
	else
	echo "$fileName $@ completed, See attachment" > mytempfile; cat mytempfile | mutt -s "$fileName $@ completed" -a email_log.html -- $email_list
fi
rm mytempfile
rm email_log.html
echo ""
echo "****************************** END ******************************"
echo ""
#
#  Created by Stuber, Tod P - APHIS on 5/3/2014.
#2015-04-20#
