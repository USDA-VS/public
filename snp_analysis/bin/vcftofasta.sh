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
echo "Start Time: `date`" > sectiontime
starttime=`date +%s`
argUsed="$1"
uniqdate=`date "+%Y-%m-%dat%Hh%Mm%Ss"`
echo "start time: $uniqdate"

####################################################
filterdir="/home/shared/${uniqdate}-FilterFiles"
mkdir ${filterdir}
FilterDirectory=${filterdir} #Files containing positions to filter

####################################################
function filterFileCreations () {
# Create "here-document"
cat >./FilterFileCreations.sh <<EOL

#!/bin/sh

# Use to make filter files from the text pasted from the Excel worksheet.
# working directory does not need to be set.
#################################################################################
#   Set variables:

# Path to txt file containing paste from Excel worksheet.
filterFile="${filterdir}/filterFile.txt"

# Number of columns in Excel worksheet
columns=`head $filterFile | awk 'BEGIN{ FS="\t"; OFS="\t" }  END {print NF}'`

# Location filter files are output to.
output="${filterdir}"

#################################################################################


let columns=columns+1
rm ${output}*
echo "Number of columns: $columns"

count=1
while [ $count -lt ${columns} ]; do
    echo ${count}
    filename=`awk -v x=$count 'BEGIN{FS=OFS="\t"}{print $x}' $filterFile | head -n1`
    echo "Filename: $filename"
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

#
#  Created by Stuber, Tod P - APHIS on 12/16/2013.
#


EOL

chmod 755 ./FilterFileCreations.sh

./FilterFileCreations.sh

rm ./FilterFileCreations.sh

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
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/abortus1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/abortus1/vcfs"
    echo "vcftofasta.sh ran as Brucella abortus bv 1, 2 or 4"
    echo "Script vcftofasta.sh ran using Brucella abortus bv 1, 2 or 4 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == mel ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Melitensis/script_dependents/Mel_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Melitensis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Melitensis/vcfs"
    echo "vcftofasta.sh ran as B. melitensis"
    echo "Script vcftofasta.sh ran using B. melitensis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == suis1 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Suis1/script_dependents/Suis1_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Suis1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Suis1/vcfs"
    echo "vcftofasta.sh ran as B. suis bv1"
    echo "Script vcftofasta.sh ran using B. suis bv1 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == suis2 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Suis2/script_dependents/suis2_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Suis2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Suis2/vcfs/"
    echo "vcftofasta.sh ran as B. suis bv2"
    echo "Script vcftofasta.sh ran using B. suis bv2 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == suis3 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Suis3/script_dependents/Suis3_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Suis3/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Suis3/vcfs"
    echo "vcftofasta.sh ran as B. suis bv3"
    echo "Script vcftofasta.sh ran using B. suis bv3 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == suis4 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Suis4/script_dependents/Suis4_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/_Brucela/Suis4/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Suis4/_VCF"
    echo "vcftofasta.sh ran as B. suis bv4"
    echo "Script vcftofasta.sh ran using B. suis bv4 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

elif [[ $1 == canis ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Canis/script_dependents/Canis_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Canis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Canis/_VCF"
    echo "vcftofasta.sh ran as B. canis"
    echo "Script vcftofasta.sh ran using B. canis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"


elif [[ $1 == ceti1 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/ceti1/script_dependents/Ceti1_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/ceti1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/ceti1/_VCF"
    echo "vcftofasta.sh ran as B ceti group 1"
    echo "Script vcftofasta.sh ran using B ceti group 1 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"


elif [[ $1 == ceti2 ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/ceti2/script_dependents/Ceti2_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/ceti2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/ceti2/_VCF"
    echo "vcftofasta.sh ran as B ceti group 2"
    echo "Script vcftofasta.sh ran using B ceti group 2 variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"


elif [[ $1 == ovis ]]; then

    genotypingcodes="/bioinfo11/TStuber/Results/brucella/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/bioinfo11/TStuber/Results/brucella/Ovis/script_dependents/Ovis_Defining_SNPs.txt"
    coverageFiles="/bioinfo11/TStuber/Results/brucella/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    RemoveFromAnalysis="/bioinfo11/TStuber/Results/brucella/Ovis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/bioinfo11/TStuber/Results/brucella/Ovis/vcfs"
    echo "vcftofasta.sh ran as B. ovis"
    echo "Script vcftofasta.sh ran using B. ovis variables" > section5
    email_list="tod.p.stuber@usda.gov Christine.R.Quance@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    if [[ $2 == elite ]]; then
        echo "Only the "elite" bovis isolates are being ran"
        sleep 5
    else
        echo "All bovis are being ran"
        echo "Like to run selected isolates? Use... vcftofasta.sh bovis elite"
        sleep 5
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb1/tb1Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb2/tb2Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb3/tb3Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4a/tb4aFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb4b/tb4bFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb5/tb5Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/tbc/tb6/tb6Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > /home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt
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
    email_list="tod.p.stuber@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"

    excelinfile="/bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/vcfs/Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' | tr -d "'"  > ${filterdir}/filterFile.txt
    filterFileCreations

else

    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, bovis, tb1, tb2, tb3, tb4a, tb4b, tb5, tb6, para"
    echo "For example, type ~$ vcftofasta.sh bovis"
    echo ""
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

Ncov=1 # Coverage below this value will be changed to -

fulDir=$PWD # Current working directory, do not change.

#################################################################################

# Count the number of chromosomes used in the reference when VCFs were made.
singleFile=`ls *.vcf | head -1`
chromCount=`awk ' $0 !~ /^#/ {print $1}' $singleFile | sort | uniq -d | awk 'END {print NR}'`
chroms=`awk ' $0 !~ /^#/ {print $1}' $singleFile | sort -n | uniq -d`
echo "The number of chromosomes seen in VCF: $chromCount"

#################################################################################

# Remove selected files from comparison
# Use file:  /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysis.txt

function removeIsolates () {

cat ${RemoveFromAnalysis} | tr '\r' '\n' | awk '{print $1}' > /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysisUnixReady.txt

removeList=`cat /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/RemoveFromAnalysisUnixReady.txt`

for i in $removeList; do
    rm *${i}*
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
        	echo "$i has data"
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

echo "AConeCallPosition.sh is running"
#echo "*********************************************************************" >> section2
#echo "Possible Mixed Isolates" > section2
#echo "Defining SNPs that are called as AC=1" >> section2
echo "" >> section2
for i in *.vcf; do
(echo "Finding possible AC1 matches for $i"; for pos in $positionList; do awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ "^"x"$" ) print FILENAME, "Pos:", $2, "QUAL:", $6, $8 }' $i; done | grep "AC=1;A" | awk 'BEGIN {FS=";"} {print $1, $2}' >> section2) &
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
        echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
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

# This function make sure there is a matching coverage file for each vcf.

function checkMatchingCoverageFile () {

echo "Checking for coverageFile for each vcf"

# Make list of tb numbers of vcfs.
for i in *; do
getbase=`basename "$i"`
number=`echo $getbase | sed $tbNumberV | sed $tbNumberW`
echo $number >> vcfFiles
done

# Make list of tb numbers of coverageFiles
for f in $coverageFiles/*; do
getbase=`basename "$f"`
number=`echo $getbase | sed $tbNumberV | sed $tbNumberW | sed 's/-coverage//g'`
echo $number >> coverageFiles
done

uniqlist=`cat vcfFiles coverageFiles | sort | uniq -u`

rm vcfFiles
rm coverageFiles

uniqNumberSize=`echo $uniqlist | wc | awk '{print $3}'`
if [ $uniqNumberSize -gt 4 ]
then
echo "No matching file was present"
echo "Either there is no coverageFile for a vcf"
echo "Or there is an extra coverageFile that is not needed"
echo "Must have matching files, and restart script."
echo "$uniqlist"

else
echo "coverageFiles are available for all vcfs"
echo "Perfect match, Good Job Chap!"
fi

}

#################################################################################

# Change SNPs with low QUAL values to N, based on parameter set above in variable settings

function changeLowCalls () {

for i in *.vcf; do
(echo "Changeing low calls"; base=`basename $i .vcf`; awk -v x=$lowEnd -v y=$highEnd 'BEGIN {OFS="\t"} { if ($6 >= x && $6 <= y) print $1, $2, $3, $4, "N", $6, $7, $8; else print $0 }' $i > ${base}.txt; rm $i; mv ${base}.txt ${base}.vcf) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

}

#################################################################################

#   Function: fasta and table creation
function fasta_table () {

# Loop through the directories
directories=`ls`
echo "$directories"

for d in $directories; do
    cd ./$d/
    dir=`basename $PWD`
    echo "Directory:  $dir"
    
	mkdir starting_files
    cp *.vcf ./starting_files
    echo "***Grabbing vcf file names"

    if [ $FilterGroups == yes ]; then
         if [ $((chromCount)) -eq 1 ]; then
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
            for i in *.vcf; do
                (m=`basename "$i"`; n=`echo $m | sed $dropEXT`
                echo "***Adding filter to $n***"
                awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' $i > $i.file
                cat "${FilterDirectory}/$d.txt" $i.file >> $i.catFile
                cat $i.catFile | sort | uniq -d > $i.txt
                pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
                echo $pos

                awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.subfilter.vcf
                
		rm $i.file
                rm $i.catFile
                rm $i.txt
                grep -v "Not_Included" $n.subfilter.vcf > $n.nomisfits.vcf
                mv $n.nomisfits.vcf $i
		)  &
                let count+=1
                [[ $((count%NR_CPUS)) -eq 0 ]] && wait
            done
        elif [ $((chromCount)) -eq 2 ]; then
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                for i in *.vcf; do m=`basename "$i"`; n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
                #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                echo "***Adding filter to $n***"
                awk 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ {print $1, $2}' $i > $i.file
                cat "${FilterDirectory}/$d.txt" $i.file >> $i.catFile
                cat $i.catFile | sort -k1,1 -k2,2 | uniq -d > $i.txt
                pos1=`cat $i.txt | grep "chrom1" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
                echo "pos1: $pos1"
                awk -v var1="chrom1" -v var2=$pos1 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.subfilter1.vcf
                pos2=`cat $i.txt | grep "chrom2" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
                echo "pos2: $pos2"
                awk -v var1="chrom2" -v var2=$pos2 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $n.subfilter1.vcf > $n.subfilter.vcf
                rm $i.file
                rm $i.catFile
                rm $i.txt
                rm $n.subfilter1.vcf
                grep -v "Not_Included" $n.subfilter.vcf > $n.nomisfits.vcf
                mv $n.nomisfits.vcf $i

                done
        else
        echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
        fi
    wait
    sleep 2

    mkdir marked_files
    mv *.subfilter.vcf ./marked_files
    else
        echo "********** Group/Subgroup/Clade Filter not ran ***********"
        echo "********** Group/Subgroup/Clade Filter not ran ***********" >> ../../section5
    fi

    # Make concatemer with the position and REF call.
    echo "***Making Concatemer"
    for i in *.vcf; do
    awk -v Q="$QUAL" '$0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' $i >> concatemer
    done

    # Get rid of duplicates in concatemer and list all the positions and REF calls
    echo "***Making total_pos"
    cat concatemer | sort -nk1 | uniq | sort -k1.6n -k1.8n > total_pos

# Find AC1 positions also found in total_pos    
    awk '{print $1}' total_pos > total.list

    for i in *.vcf; do 
	(m=`basename "$i"`; n=`echo $m | sed 's/\..*//'`
	# search for AC1 positions
	awk ' $0 !~ /^#/ && $8 ~ /^AC=1/ && $6 > 0 {print $1 "-" $2}' $i > ${n}.list
	# AC1 positions that are being found in this group
	positionsfound=`cat ${n}.list total.list | sort -n | uniq -d`
	countfind=`echo $positionsfound | wc -w`
	echo "positonsfound: $positionsfound  countfind: $countfind"
    	rm ${n}.list
	if [[ -z $positionsfound ]]; then
		positionsfound="No positions found"
	fi
	
	if [[ $countfind -gt 2  ]]; then
	searchname=`echo $n | sed 's/_.*//'`

        	if [[  $argUsed == para ]]; then
            	unmappedContigs=`grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/data/${searchname}/BWAmem-GATK/QualityValues/*stats.txt`
        	elif [[  $argUsed == bovis ]]; then
            	unmappedContigs=`grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script1/${searchname}/BWAmem-GATK/QualityValues/*stats.txt`
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

    totalSNPs=`grep -c ".*" total_pos`
    echo "$d total SNPs: $totalSNPs" >> ../../section4

#    echo "***The total_pos of $d" >> ../../section4
#    grep -c ".*" total_pos >> ../../section4

    echo "***Creating normalized vcf using AC2, QUAL > $QUAL"
        # Grab the name of the vcf file

#########################################################################
for i in *.vcf; do
    (n=${i%.vcf}
    # echo the name grabbed
    echo $n
    # Create .cut file that lists the positions and ALT calls
    awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' $i > $n.cut

    # Fill in the .cut file with REF calls at positions that were not called as SNPs
    #cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1 > $n.filledcutnoN
    cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -k1.6n -k1.8n > $n.filledcutnoN

    ##############################################################
    # Change zero coverage regions

    # get positions being used
    awk '{print $1}' $n.filledcutnoN > ${n}.filledcutNumbers

    # Get zero coverage positions.
    awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' ${i} > ${n}.zeropositions

    # zero duplicate positions will need to be kept
    cat ${n}.filledcutNumbers ${n}.zeropositions | sort | uniq -d > ${n}.zerotokeep

    # get zero position, these are only positions already in the filledcutnoN
    if [ -s ${n}.zerotokeep ]; then
        grep -f ${n}.zerotokeep ${n}.zeropositions | awk '{print $1, "-"}'> ${n}.zerotomerge
        # merge zero updates to filledcut
        cat ${n}.zerotomerge $n.filledcutnoN | awk '{ if (a[$1]++ == 0) print $0; }' | sort -k1.6n -k1.8n > ${n}.filledcut
        rm ${n}.filledcutnoN
        rm ${n}.zerotomerge
    else
        mv $n.filledcutnoN ${n}.filledcut
    fi

    rm ${n}.filledcutNumbers
    rm ${n}.zeropositions
    rm ${n}.zerotokeep)  &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

#########################################################################

echo "sleeping 5 seconds at line number: $LINENO"; sleep 5
wait

        # Make a concatemer of the .filledcut files
        for i in *.filledcut; do
            cat $i >> cutConcatemer
            done
        echo "***Making the select file containing positions of interest"
        # Capture only positions that have more than one SNP type called at a position
        cat cutConcatemer | sort -k1.6n -k1.8n | uniq | awk '{print $1}' | uniq -d > select
        # Compare the positions in select with total_pos and output total_pos position that match select positions
        # but only with positions that are in the select file.
        # This getting rid of calls that are the same for all isolates being analyzed
        echo "***grepping the total_pos file"

        echo "***grepping the .filledcut files for $d"

        grep -f select total_pos | sort -k1.6n -k1.8n > clean_total_pos
        
        # Begin the table
        awk '{print $1}' clean_total_pos | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> $d.table.txt
        awk '{print $2}' clean_total_pos | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}'>> $d.table.txt
                    # Make the fasta files:  Fill in positions with REF if not present in .clean file

        for i in *.filledcut; do
            (m=`basename "$i"`
            n=`echo $m | sed $dropEXT`
            # Compare the positions in select with "isolate".cut and output position for .cut that only matched select positions
            egrep -f select $i | sort -k1.6n -k1.8n > $n.pretod

        ##############################################################
        # Change AC1s to IUPAC

            # get positions being used
            awk '{print $1}' ${n}.pretod > ${n}.usedpostions
            # get AC1 positions and iupac calls  that were changed to iupac
            awk ' $0 !~ /#/ && $6 > 300 && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' ${i%filledcut}vcf > ${n}.ac
            # get just positions of those AC1 grabbed above
            awk '{print $1}' ${n}.ac > ${n}.acpositions
           # AC duplicate positions will need to be kept
            cat ${n}.usedpostions ${n}.acpositions | sort | uniq -d > ${n}.actokeep
            # get AC1 position with iupac, these are only positions already in the pretod

            if [ -s ${n}.actokeep ]; then
                grep -f ${n}.actokeep ${n}.ac > ${n}.actomerge
                # merge iupac updates to filledcut
                cat ${n}.actomerge $n.pretod | awk '{ if (a[$1]++ == 0) print $0; }' | sort -k1.6n -k1.8n > $n.tod
                rm ${n}.pretod
                rm ${n}.actomerge
            else
                mv $n.pretod $n.tod
            fi
            rm ${n}.usedpostions
            rm ${n}.ac
            rm ${n}.acpositions
            rm ${n}.actokeep
            rm ${i%filledcut}vcf
        ##############################################################

	sed 's/chrom[0-9-]*//g' $n.tod | tr -d [:space:] | awk '{print $0}' | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g'  > $n.fas
            # Add each isolate to the table
            awk '{print $2}' $n.tod | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> $d.table.txt) &
    		let count+=1
    		[[ $((count%NR_CPUS)) -eq 0 ]] && wait
	done
wait

echo "sleeping 5 seconds at line number: $LINENO"; sleep 5
        #Make a reference fasta sequence
        awk '{print $2}' clean_total_pos > root
        cat root | tr -cd "[:print:]" | sed "s/^/>root;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > root.fas
	echo "" >> root.fas

	totalSNPs=`grep -c ".*" total_pos`
	echo "$d informative SNPs: $totalSNPs" >> ../../section4

        # Make a file containing all fasta files.
        cat *.fas > ${dir}_alignment.fasta

        #Clean-up
        echo "***Cleaning folder"
        rm *.cut
        rm *.filledcut
        rm concatemer
        rm cutConcatemer
        rm *.tod
        mkdir fasta
        mv *.fas ./fasta
        rm total_pos
        rm select
        rm root
        rm clean_total_pos
        cp /home/shared/Table_Template.xlsx ./${d}-Table_Template.xlsx
	echo "***Done"
    cd ..
    done
echo "sleeping 5 seconds at line number: $LINENO"; sleep 5

}

#****************************************************************
function alignTable () {

# Beginning in fasta folder
echo "$d *********"
pwd

cat *.fas | sed '/root/{N;d;}' >> fastaGroup.txt
cat *.fas >> RAxMLfastaGroup.txt

#clustalw2 -OUTFILE=alignment.txt -RANGE=1,2 -OUTPUT=FASTA -INFILE=fastaGroup.txt & 
/usr/local/bin/standard-RAxML-master/raxmlHPC-SSE3 -s RAxMLfastaGroup.txt -n ${d} -m GTRCAT -p 12345 && nw_reroot RAxML_bestTree.${d} root | nw_display -s -w 1000 -v 20 -b 'opacity:0' -i 'font-size:8' -l 'font-family:serif;font-style:italic' -d 'stroke-width:2;stroke:blue' - > ../${d}-tree.svg && inkscape -f ../${d}-tree.svg -A ../${d}-tree.pdf &
wait
rm RAxML_parsimonyTree*
for i in RAxML*Tree*; do mv $i ../${i}.tre; done
#grep ">" alignment.txt | sed 's/>//g' > cleanedAlignment.txt

pwd
inputfile=`ls RAxML_result*`
tr ":" "\n" < ${inputfile} | tr "," "\n" | sed 's/(//g' | sed 's/)//g' | grep -v "\.[0-9]*" | grep -v "root" > cleanedAlignment.txt

awk 'NR==FNR{o[FNR]=$1; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' cleanedAlignment.txt ../$d.table.txt > joined.txt
grep "reference" ../$d.table.txt > references.txt
cat references.txt joined.txt >> joined2.txt
mv joined2.txt ../$d.sortedTable.txt

#rm alignment.txt
#rm cleanedAlignment.txt
#rm *.dnd
#rm fastaGroup.txt
rm joined.txt
rm references.txt

echo "**** orgTable.sh Started ****"
cd ..
#get the number of columns
columnCount=`awk '$0 !~ /^$/ {print $0}' *sortedTable.txt | awk '{print NF-1; exit}'`
echo "Column Count: $columnCount"

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
#awk -f /Users/Shared/_programs/_my_scripts/awk_scripts/transpose.awk orginizedTable6.txt > orginizedTable7.txt
#awk '{for (i=1; i<=NF; i++)  { a[NR,i] = $i}} NF>p { p = NF } END {for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]} print str}} orginizedTable6.txt > orginizedTable7.txt
###
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
###

#Orgainize file based on 1st 2 columns
sort -n -k1 orginizedTable7.txt | sort -n -k2 > orginizedTable8.txt

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable8.txt > orginizedTable9.txt
awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' orginizedTable9.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > orginizedTable10.txt

#Transpose back
awk -f /home/tstuber/workspace/stuber/awk_scripts/transpose.awk orginizedTable10.txt > orginizedTable11.txt

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
echo "**** orgTable.sh Finished ****"

}
#****************************************************************


################################################################################
#################################################################################
#################################################################################
echo "****************************** START ******************************"
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

#Prepare Filter files.
filterFilespreparation

#Test for match coverage file
#checkMatchingCoverageFile

#copy the original vcfs to /starting_files
mkdir starting_files
mv *.* ./starting_files

# If bovis are ran default will only run with files check "misc" in FileMaker
# Untitled.tab exported from FileMaker must contain "isolate names" followed by "Misc".

	if [[ $2 == elite ]]; then
        echo "Only analyzing elite files"

        for i in `cat elite`; do
        name=`ls starting_files | grep $i`
        cp ./starting_files/$name ./
        done

        for i in `find ./starting_files/ -mtime -30`; do
        cp $i ./
        done

    else
        echo "all samples will be ran"
        cp ./starting_files/* ./
	fi
rm elite

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
    echo "VALtest: $VALtest"
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
            echo "A"
        else
            name=${searchName}-Val
            mv $i ${name}.vcf
            echo "B"
        fi
    else # New name WAS found
        if [[ -z "$VALtest" ]]; then
            name=$newName
            mv $i ${name}.vcf
            echo "C"
        else
            name=${newName}-Val
            echo "newName is $name"
            mv $i ${name}.vcf
            echo "D"
        fi
    fi
done
rm outfile
##################### Start: Make Files Unix Compatiable #####################

#Fix validated (VAL) vcf files.  This is used in vcftofasta scripts to prepare validated vcf files opened and saved in Excel.
#Create list of isolates containing "VAL"
#Do NOT make this a child process.  It messes changing column 1 to chrom
echo "Making Files Unix Compatiable"
for v in *.vcf; do
    (dos2unix $v #Fixes files opened and saved in Excel
    cat $v | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' > $v.temp
    mv $v.temp $v) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

############## Capture the number of chromosomes and their name from a single VCF ##############

echo "The chromosome count is: $chromCount"
pwd

if [ $chromCount -gt 2 ]; then
    echo "More than 2 chromosomes detected."
    echo "Unable to handle more than 2."
    echo "Exiting script at line $LINENO"
    exit 1
fi
# Change chromosome identification to general chrom1 and/or chrom2

for f in *.vcf; do
    echo "echoing f: $f"
    num=1
    for i in $chroms; do
        sed "s/$i/chrom${num}/g" $f > temp.vcf
        mv temp.vcf $f
        echo "$i was marked as chrom${num}"
    num=$(( $num + 1 ))
    done
done

########################################################################

printf "%s\t%s\t%s\t%s\n" "TB Number" "Group" "Subgroup" "Clade" > FileMakerGroupImport.txt

AConeCallPosition
wait

# Change low QUAL SNPs to N, see set variables
changeLowCalls
wait

######################## Change AC1s to IUPAC ########################

for i in *vcf; do

(echo "Changing AC=1 to IUPAC: $i"
awk '
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
echo "***Marking all VCFs and removing filtering regions"
#echo "***Marking all VCFs and removing filtering regions was done." >> log
    # Label filter field for positions to be filtered in all VCFs
        if [ $((chromCount)) -eq 1 ]; then
        for i in *.vcf; do
        (m=`basename "$i"`; n=`echo $m | sed $dropEXT`; echo "********* $n **********"
        # Get usable positions in the VCF
        awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' $i > $i.file
        # Combine with positions that will be filtered
        cat "${FilterDirectory}/FilterToAll.txt" $i.file >> $i.catFile
        
# Output any duplicate positions, aka decreasing positions to be marked and used by awk
        cat $i.catFile | sort | uniq -d > $i.txt
        # preparing postions
        pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`

	# Making a vcf with positions marked that should not be included based on filter file
        awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filter.vcf
        # Clean up
	
        rm $i.file; rm $i.catFile; rm $i.txt
        # Removed positions
        grep -v "Not_Included" $n.filter.vcf > $n.noPPE.vcf
        # Change name
        mv $n.noPPE.vcf $i) &
     	let count+=1
      	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
        done

        elif [ $((chromCount)) -eq 2 ]; then
        for i in *.vcf; do

            m=`basename "$i"`
            n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
            echo "********* $n **********"
            awk 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ {print $1, $2}' $i > $i.file
            cat "${FilterDirectory}/FilterToAll.txt" $i.file >> $i.catFile
cat $i.catFile | sort -k1,1 -k2,2 | uniq -d > $i.txt
            pos1=`cat $i.txt | grep "chrom1" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
            echo "pos1: $pos1"
            awk -v var1="chrom1" -v var2=$pos1 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filter1.vcf
            pos2=`cat $i.txt | grep "chrom2" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
            echo "pos2: $pos2"
            awk -v var1="chrom2" -v var2=$pos2 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $n.filter1.vcf > $n.filter.vcf
           
	rm $i.file
            rm $i.catFile
            rm $i.txt
            rm $n.filter1.vcf
            grep -v "Not_Included" $n.filter.vcf > $n.noPPE.vcf
            mv $n.noPPE.vcf $i
        done
        else
            echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
        fi

    wait
    sleep 2

        mkdir marked_files
        mv *.filter.vcf ./marked_files
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

    # Get quality positions in VCF
    formatedpos=`awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/{print $2}' $i | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/' | sed 's/|$$//'`
    
    # If a group number matches a quality position in the VCF (formatedpos) then print the position
    groupNumbers=`grep "Group" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Group Numbers: $groupNumbers"

    # Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
    sizeGroup=`echo $groupNumbers | wc | awk '{print $3}'`

    # Loop through the number of groups positions found
    loops=`grep "Group" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`
    if [ $sizeGroup -lt 2 ]; then # There was not a position found that places VCF into group
        echo "$i Grp not found" >> section3
        echo "$i was not assigned a Group"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: all_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to all_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                        echo "$i Multi Grp" >> section3
                    fi
                echo "The counter is at $counter"
                groupNumber=`echo "$groupNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $groupNumber
                mkdir -p Group-$groupNumber #Make groupNumber folder if one does not exist.
                cp $i ./Group-$groupNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

    # Get the Subgroup number
    subgroupNumbers=`grep "Subgroup" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Subgroup Numbers: $subgroupNumbers"
        # Check if a Group, subgroup or clade was called.
        sizeGroup=`echo $subgroupNumbers | wc | awk '{print $3}'`

    loops=`grep "Subgroup" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`

        if [ $sizeGroup -lt 2 ];
        then
        echo "$i was not assigned a Subgroup"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: all_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to all_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                    echo "$i Multi Subgrps" >> section3
                    fi
                echo "The counter is at $counter"
                subgroupNumber=`echo "$subgroupNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $subgroupNumber
                mkdir -p Subgroup-$subgroupNumber #Make groupNumber folder if one does not exist.
                cp $i ./Subgroup-$subgroupNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

    # Get the Clade number
    cladeNumbers=`grep "Clade" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Clade Numbers: $cladeNumbers"

echo "${i%.vcf} $groupNumbers $subgroupNumbers $cladeNumbers" >> section3
tbn=`echo $i | sed 's/_.*//'`
printf "%s\t%s\t%s\t%s\n" "${tbn}" "$groupNumbers" "$subgroupNumbers" "$cladeNumbers" >> FileMakerGroupImport.txt
        # Check if a Group, subgroup or clade was called.
        sizeGroup=`echo $cladeNumbers | wc | awk '{print $3}'`
        loops=`grep "Clade" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`
        if [ $sizeGroup -lt 2 ];
        then

        echo "$i was not assigned a Clade"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: all_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to all_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                    echo "$i Multi Clade" >> section3
                    fi
                echo "The counter is at $counter"
                cladeNumber=`echo "$cladeNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $cladeNumber
                mkdir -p Clade-$cladeNumber #Make groupNumber folder if one does not exist.
                cp $i ./Clade-$cladeNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

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
# Make concatemer with the position and REF call.
    # Factor in possible multiple chromosomes
    echo "***Making Concatemer"
    for i in *.vcf; do
awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' $i >> concatemer
    done

# Get rid of duplicates in concatemer and list all the positions and REF calls
echo "***Making total_pos"
cat concatemer | sort -nk1 | uniq | sort -k1.6n -k1.8n > total_pos

# Count the number of SNPs
totalSNPs=`grep -c ".*" total_pos`
echo "Total SNPs: $totalSNPs" >> ../section4
#grep -c ".*" total_pos >> ../section4
echo "***Creating normalized vcf using AC2, QUAL > 150"
# Grab the name of the vcf file

for i in *.vcf; do
    (n=${i%.vcf}
    echo $n
    awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' > $n.cut

    #cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1 > $n.filledcutnoN
    cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -k1.6n -k1.8n > $n.filledcutnoN

    ##############################################################
    # Change zero coverage regions
    echo "Change zero coverage regions"
    # get positions being used
    awk '{print $1}' $n.filledcutnoN > ${n}.filledcutNumbers

    # Get zero coverage positions.
    awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' ${i} > ${n}.zeropositions

    # zero duplicate positions will need to be kept
    cat ${n}.filledcutNumbers ${n}.zeropositions | sort | uniq -d > ${n}.zerotokeep

    # get zero position, these are only positions already in the filledcutnoN
    if [ -s ${n}.zerotokeep ]; then
    grep -f ${n}.zerotokeep ${n}.zeropositions | awk '{print $1, "-"}' > ${n}.zerotomerge
    # merge zero updates to filledcut
    cat ${n}.zerotomerge $n.filledcutnoN | awk '{ if (a[$1]++ == 0) print $0; }' | sort -k1.6n -k1.8n > ${n}.filledcut
    rm ${n}.filledcutnoN
    rm ${n}.zerotomerge
    else
    mv $n.filledcutnoN ${n}.filledcut
    fi
    
    rm ${n}.filledcutNumbers
    rm ${n}.zeropositions
    rm ${n}.zerotokeep)  &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
wait

#########################################################################

echo "sleeping 5 seconds at line number: $LINENO"; sleep 5
wait

# Begin the table
awk '{print $1}' total_pos | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> all_vcfs.table.txt
awk '{print $2}' total_pos | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> all_vcfs.table.txt
echo "***grepping the .filledcut files"
# Make the fasta files:  Fill in positions with REF if not present in .clean file

for i in *.filledcut; do
    (echo " working on filled cut for $i"
    m=`basename "$i"`
    n=`echo $m | sed $dropEXT`
    # Compare the positions in select with "isolate".cut and output position for .cut that only matched select positions
    #grep -f select $i > $n.tod

    # Use this cat command to skip the time intensive grep
    # With all_vcfs this grep doesn't eliminate many snps.
    sed 's/chrom[0-9-]*//g' $i | tr -d [:space:] | awk '{print $0}' | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g'  > $n.fas

    # Add each isolate to the table
    awk '{print $2}' $i | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> all_vcfs.table.txt) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
wait
wait

done
wait
echo "sleeping 5 seconds at line number: $LINENO"; sleep 5

echo "grepping filledcut files is finished"
#Make a reference fasta sequence
awk '{print $2}' total_pos > root
cat root | tr -cd "[:print:]" | sed "s/^/>root;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > root.fas
echo "" >> root.fas

totalSNPs=`grep -c ".*" total_pos`
echo "Total informative SNPs: $totalSNPs" >> ../section4

#Clean-up
mkdir starting_files
echo "***Cleaning folder"
mv *.vcf ./starting_files
rm *.cut
rm *.filledcut
rm *.filledcutnoN
rm concatemer
rm cutConcatemer
rm *.tod
mkdir fasta
mv *.fas ./fasta
#rm total_pos
rm root

if [[ $2 == elite ]]; then
	d="all_vcfs"
        cd ./fasta
        alignTable
else
	echo "Tree not made when all samples are ran"	
fi

echo "***Done"
echo "Full Directory: ${fulDir}"

##################### End: All vcf folder #####################

echo "***************************************************"
echo "***************** STARTING Groups *****************"
echo "***************************************************"
# Change directory to all_groups
cd ${fulDir}/all_groups
fasta_table &

echo "***************************************************"
echo "**************** STARTING SUBGROUPS ***************"
echo "***************************************************"
# Change directory to all_subgroups
cd ${fulDir}/all_subgroups
fasta_table &

echo "***************************************************"
echo "***************** STARTING CLADES *****************"
echo "***************************************************"
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

echo "***************************************************"
echo "********** STARTING all_clades Alignment **********"
echo "***************************************************"
cd ${fulDir}/all_clades

workingdir=`basename $PWD`

if [ $workingdir == all_clades ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_clades/${d}/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
	alignTable & 

    pwd
done

else
echo "*** $workingdir not found ***"
fi

echo "***************************************************"
echo "********** STARTING all_groups Alignment **********"
echo "***************************************************"
cd ${fulDir}/all_groups

workingdir=`basename $PWD`

if [ $workingdir == all_groups ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_groups/${d}/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
    alignTable &
pwd
done
else
echo "*** $workingdir not found ***"
fi

echo "***************************************************"
echo "******** STARTING all_subgroups Alignment *********"
echo "***************************************************"
cd ${fulDir}/all_subgroups

workingdir=`basename $PWD`

if [ $workingdir == all_subgroups ]
then

directories=`ls`
for d in $directories; do
    cd ${fulDir}/all_subgroups/$d/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
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
rm -r ${FilterDirectory}

echo "Copying to ${bioinfoVCF}"
cp -r $PWD ${bioinfoVCF}
echo "******* $LINENO, $PWD"
fileName=`basename $0`

if [[ $3 == me ]]; then
    email_list="Tod.P.Stuber@aphis.usda.gov"
	echo "vcftofasta.sh completed" > mytempfile; cat mytempfile | mutt -s "vcftofasta.sh completed subject" -a email_log.html -- $email_list
	else
	echo "$fileName $@ completed, See attachment" > mytempfile; cat mytempfile | mutt -s "$fileName $@ completed" -a email_log.html -- $email_list
fi
rm mytempfile
rm email_log.html

echo "****************************** END ******************************"
pwd

#
#  Created by Stuber, Tod P - APHIS on 5/3/2014.
#2015-04-20#
