#!/bin/sh

#  bovis_processzips.sh
#  Working directory should contain paired-end fastq reads
#  Reads must be included as _R1 and _R2
#  See loopfiles.sh and email_loopfiles for multiple samples.

#################################################################################
#  Dependencies ---
#   bwa, http://bio-bwa.sourceforge.net/bwa.shtml
#   samtools, http://$HOME/prog/samtools-1.3/samtools.sourceforge.net/$HOME/prog/samtools-1.3/samtools.shtml
#   picard, http://picard.sourceforge.net/command-line-overview.shtml
#   gatk, http://www.broadinstitute.org/gatk/
#   bamtools
#   File containing high quality SNPs, Volumes/Mycobacterium/Go_To_File/HighestQualitySNPs.vcf
#   Reference in fasta format, /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta
#################################################################################

cpu=$(($(nproc) - 1))

echo "**************************************************************************"
echo "**************************** START ${PWD##*/} ****************************"
echo "**************************************************************************"

#debug stuff
alias pause='read -p "$LINENO Enter"'

picard="$HOME/prog/picard-tools-2.1.0/picard.jar"
gatk="$HOME/prog/gatk/GenomeAnalysisTK.jar"
igvtools="$HOME/prog/IGVTools/igvtools.jar"
BRUC_MLST=$(which Bruc_MLST.sh)
SPOLIGOSPACERFINDER=$(which spoligoSpacerFinder.sh)

echo "current directory"
pwd
startingdir=$(pwd)

#script dependenties
dependents="$HOME/pipeline/public-master/snp_analysis/script_dependents"

# Move zip files to their own directory
mkdir ./zips
mv *.fastq* ./zips
mkdir ./bwamem-gatk
cd bwamem-gatk/

# Make alias links in BWA-GATK directory to zip files
ls ../zips/*.fastq* | while read file; do ln -s $file; done

if [ $1 == ab1 ]; then
    cp $dependents/Brucella_abortus/NC_00693c.fasta ./
    hqs="$dependents/Brucella_abortus/NC_00693cHighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == mel ]; then
    cp $dependents/Brucella_melitensis/BmelitensisM5-90.fasta ./
    hqs="$dependents/Brucella_melitensis/melHighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == suis1 ]; then
    cp $dependents/Brucella_suis_bv1/NC_01725c.fasta ./
    hqs="$dependents/Brucella_suis_bv1/NC_01725cHighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == suis2 ]; then
    cp $dependents/Brucella_suis_bv2/Bsuisbv2-94-11.fasta ./
    hqs="$dependents/Brucella_suis_bv2/suis2HighestQualitySNPs.vcf"
    
    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == suis3 ]; then
    cp $dependents/Brucella_suis_bv3/B-REF-BS3-686.fasta ./
    hqs="$dependents/Brucella_suis_bv3/suis3HighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == suis4 ]; then
    cp $dependents/Brucella_suis_bv4/B-REF-BS4-40.fasta ./
    hqs="$dependents/Brucella_suis_bv4/suis4HighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == canis ]; then
    cp $dependents/Brucella_canis/BcanisATCC23365.fasta ./
    hqs="$dependents/Brucella_canis/canisHighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == ceti1 ]; then
    cp $dependents/Brucella_ceti-grp1/Bceti1Cudo.fasta ./
    hqs="$dependents/Brucella_ceti-grp1/ceti1HighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == ceti2 ]; then
    cp $dependents/Brucella_ceti-grp2/Bceti2-TE10759.fasta ./
    hqs="$dependents/Brucella_ceti-grp2/ceti2HighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $1 == ovis ]; then
    cp $dependents/Brucella_ovis/BovisATCC25840.fasta ./
    hqs="$dependents/Brucella_ovis/BovisATCC25840HighestQualitySNPs.vcf"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../bwamem-gatk/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

# Lineage 1
elif [ $1 == TB1 ]; then
    cp $dependents/TB1/NC_017528.fasta ./
    hqs="$dependents/TB1/HQ-NC_017528.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 2
elif [ $1 == TB2 ]; then
    cp $dependents/TB2/NC_021251.fasta ./
    hqs="$dependents/TB2/HQ-NC021251.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 3
elif [ $1 == TB3 ]; then
    #cp $dependents/TB3/NC_021193.fasta ./
    #hqs="$dependents/TB3/HQ-13-7575.vcf"
    cp $dependents/TB3/NC_021193it3-readreference.fasta ./
    hqs="$dependents/TB3/13-7575-highqualitysnps.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 4.1 and 4.2
elif [ $1 == TB4a ]; then
    cp $dependents/TB4a/NC002755.fasta ./
    hqs="$dependents/TB4a/HQ-NC002755.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 4.9
elif [ $1 == TB4b ]; then
    cp $dependents/TB4b/NC018143.fasta ./
    hqs="$dependents/TB4b/HQ-NC018143.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 5
elif [ $1 == TB5 ]; then
    cp $dependents/TB5/APKD01000001.fasta ./
    hqs="$dependents/TB5/HQ-16-2185-11.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage 6
elif [ $1 == TB6 ]; then
    cp $dependents/TB6/NC_015758.fasta ./
    hqs="$dependents/TB6/HQ-NC015758.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage Bov-Afri, AF2122
elif [ $1 == TBBOV ]; then
    cp $dependents/Mycobacterium_bovis/NC_002945.fasta ./
    hqs="$dependents/Mycobacterium_bovis/HighestQualitySNPs.vcf"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

elif [ $1 == para ]; then
   cp $dependents/paraTB/NC_002944.fasta ./
   hqs="$dependents/paraTB/HQ-NC002944.vcf"

    ###################################################################

else
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, TB1, TB2, TB3, TB4a, TB4b, TB5, TB6, TBBOV, para"
    exit 1
fi

#Reference genome
ref=$(ls | grep .fasta)
echo "Reference Input:  $ref"

#Forward reads
forReads=$(ls | grep _R1)
echo "Forward Reads:  $forReads"

#Reverese reads
revReads=$(ls | grep _R2)
echo "Reverse Reads:  $revReads"

#Retrieves reference name and name from sorted BAM file name
r=$(echo $ref | sed 's/\..*//')
n=$(echo $revReads | sed -e 's/_.*//' -e 's/\..*//')

echo "***Reference naming convention:  $r"
echo "***Isolate naming convention:  $n"

#ref genome indexing
echo "***Reference Genome $r is being indexed"
$HOME/prog/samtools-1.3/samtools faidx "$ref"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
$HOME/prog/bwa/bwa index "$ref"

# -t sets the number of threads/cores
# -r ST	 Specify the read group in a format like ‘@RG\tID:foo\tSM:bar’ Needed for GATK
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
echo "***Making BAM file"
rg="@RG\tID:${n}\tPL:ILLUMINA\tPU:${n}_RG1_UNIT1\tLB:${n}_LIB1\tSM:${n}"
$HOME/prog/bwa/bwa mem -M -r 1 -t "$cpu" -R "$rg" "$ref" "$forReads" "$revReads" | \
sambamba view -t "$cpu" -f bam -h -S /dev/stdin | \
sambamba sort -t "$cpu" -o ${n}.all.bam /dev/stdin # mapped and unmapped reads

#https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax
#only mapped
sambamba view -t "$cpu" -f bam -h -F "not (unmapped or mate_is_unmapped)" ${n}.all.bam | \
sambamba sort -t "$cpu" -o ${n}.sorted.bam /dev/stdin
#only unmapped
sambamba view -t "$cpu" -f bam -h -F "unmapped or mate_is_unmapped" ${n}.all.bam | \
sambamba sort -t "$cpu" -o ${n}.unmappedReads.bam /dev/stdin

#Create fastqs of unmapped reads to assemble
bedtools bamtofastq -i ${n}.unmappedReads.bam \
-fq ${n}-unmapped_R1.fastq \
-fq2 ${n}-unmapped_R2.fastq

#compress fastq
pigz ${n}-unmapped_R?.fastq

#SPAdes
# easy_install --install-dir=/home/CFIA-ACIA/duceppem/local/lib/python2.7/site-packages/ regex
spades.py -t "$cpu" -k 63 --careful -1 ${n}-unmapped_R1.fastq.gz -2 ${n}-unmapped_R2.fastq.gz -o spades_output

#file stuff
mkdir ../unmappedReads
mv ${n}-unmapped_R?.fastq.gz ../unmappedReads
mv spades_output/scaffolds.fasta ../unmappedReads
mv spades_output/spades.log ../unmappedReads
######################

echo "***Removing Duplicates and Indexing Bam"
sambamba markdup -r -t "$cpu" ${n}.sorted.bam ${n}.dup.bam
sambamba index -t "$cpu" ${n}.dup.bam

# Creates file that is used in the next step
# locally realign reads such that the number of mismatching bases is minimized across all the reads
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html
echo "***Realigner Target Creator"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T RealignerTargetCreator -I ${n}.dup.bam -R "$ref" -o ${n}.forIndelRealigner.intervals -nt "$cpu"

if [ ! -e ${n}.forIndelRealigner.intervals ]; then
	$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T RealignerTargetCreator --fix_misencoded_quality_scores -I ${n}.dup.bam -R "$ref" -o ${n}.forIndelRealigner.intervals -nt "$cpu"
fi

# Uses the RealignerTargetCreator output file to improve BAM alignment
# http://www.broadinstitute.org/gatk/guide/tagged?tag=indelrealigner
echo "***Target Intervals"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T IndelRealigner -I ${n}.dup.bam -R "$ref" -targetIntervals ${n}.forIndelRealigner.intervals -o ${n}.realignedBam.bam

if [ ! -e ${n}.realignedBam.bam ]; then
	echo "$n RealignedBam.bam failed to make.  Possible cause: Error in quality scores.  Try --fix_misencoded_quality_scores"
	echo "$n RealignedBam.bam failed to make.  Possible cause: Error in quality scores.  Try --fix_misencoded_quality_scores" > ${n}.errorReport
	#cat ${n}.errorReport | mutt -s "$n Alignment failure" -- tod.p.stuber@usda.gov
	$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T IndelRealigner --fix_misencoded_quality_scores -I ${n}.dup.bam -R "$ref" -targetIntervals ${n}.forIndelRealigner.intervals -o ${n}.realignedBam.bam
fi

# Uses a .vcf file which contains SNP calls of known high value to recalibrates base quality scores
# http://www.broadinstitute.org/gatk/guide/tagged?tag=baserecalibrator
echo "***Base Recalibrator"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T BaseRecalibrator -I ${n}.realignedBam.bam -R "$ref" -knownSites "$hqs" -o ${n}.recal_data.grp -nct "$cpu"

if [ ! -e ${n}.recal_data.grp ]; then
	$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T BaseRecalibrator --fix_misencoded_quality_scores -I ${n}.realignedBam.bam -R "$ref" -knownSites "$hqs" -o ${n}.recal_data.grp -nct "$cpu"
fi

# Make the finished "ready" .bam file
echo "***Print Reads"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T PrintReads -R "$ref" -I ${n}.realignedBam.bam -BQSR ${n}.recal_data.grp -o ${n}.preready-mem.bam -nct "$cpu"

if [ ! -e ${n}.preready-mem.bam ]; then
	$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -T PrintReads --fix_misencoded_quality_scores -R "$ref" -I ${n}.realignedBam.bam -BQSR ${n}.recal_data.grp -o ${n}.preready-mem.bam -nct "$cpu"
fi

$HOME/prog/jdk1.8.0_72/bin/./java -jar ${gatk} -T ClipReads -R "$ref" -I ${n}.preready-mem.bam -o ${n}.ready-mem.bam -filterNoBases -dcov 10
sambamba index -t "$cpu" ${n}.ready-mem.bam

#Collect Depth of coverage info
echo "***Collect Depth of Coverage"
$HOME/prog/jdk1.8.0_72/bin/./java -jar ${gatk} -T DepthOfCoverage -R "$ref" -I ${n}.preready-mem.bam -o ${n}.coverage -omitIntervals --omitLocusTable --omitPerSampleStats -nt "$cpu"

#########################

# http://www.broadinstitute.org/gatk/guide/tagged?tag=unifiedgenotyper
# In group tb4b position 3336835 was not being found in some isolates.  Adding this advance flag allowed these positions to be found.
# ploidy 2 is default

echo "***HaplotypeCaller, aka calling SNPs"
#-allowNonUniqueKmersInRef
#Can't use multithread if outputing BAM file with HaplotypeCaller!!!
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -R "$ref" -T HaplotypeCaller -I ${n}.ready-mem.bam -o ${n}.hapreadyAll.vcf -bamout ${n}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${igvtools} index ${n}.hapreadyAll.vcf

echo "******Awk VCF leaving just SNPs******"
awk '/#/ || $4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ {print $0}' ${n}.hapreadyAll.vcf > ${n}.hapreadyOnlySNPs.vcf

#Split header lines from position calls
grep "#" ${n}.hapreadyOnlySNPs.vcf > ${n}.header
grep -v "#" ${n}.hapreadyOnlySNPs.vcf > ${n}.body

#SNP positons that will be used
awk '{print $1 "%" $2}' ${n}.body > ${n}.calledSNPpositions

#Zero coverage positions
awk 'BEGIN {FS="[:\t]"} $3 == 0 {print $1 "%" $2}' ${n}.coverage > ${n}.zeroCoveragePositions
#Remove zero coverage positions that will be are in ${n}.hapreadyOnlySNPs.vcf
cat ${n}.calledSNPpositions ${n}.zeroCoveragePositions | sort | uniq -d > ${n}.duplicates
cat ${n}.zeroCoveragePositions ${n}.duplicates | sort | uniq -u > ${n}.keepTheseZeroCovPositions
zeroposition=$(grep -c ".*" ${n}.keepTheseZeroCovPositions)
refsize=$(wc -m $ref | awk '{print $1}')

#Fromat ${n}.keepTheseZeroCovPositions to VCF
sed 's/%/ /' ${n}.keepTheseZeroCovPositions | awk 'BEGIN{OFS="\t"}{print $1, $2, ".", ".", ".", ".", ".", ".", "GT", "./."}' > ${n}.vcfFormated
cat ${n}.body ${n}.vcfFormated | awk 'BEGIN{OFS="\t"}{if ($4 == ".") print $1, $2, $3, "N", $5, $6, $7, $8, $9, $10; else print $0}' > ${n}.SNPsMapzeroNoHeader.vcf
cat ${n}.header ${n}.SNPsMapzeroNoHeader.vcf > ${n}.unsortSNPsZeroCoverage.vcf
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${igvtools} sort ${n}.unsortSNPsZeroCoverage.vcf ${n}.SNPsZeroCoverage.vcf
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${igvtools} index ${n}.SNPsZeroCoverage.vcf

# Emit all sites to VCF, not just the SNPs and indels.  This allows making a UnifiedGenotyper VCF similar to what was used before using the Haplotypecaller.
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I ${n}.ready-mem.bam -o ${n}.allsites.vcf -nt "$cpu"

# This removes all positions same as the reference.  These positions are found by removing rows were column (field) 8 begins with AN=2.  This decreases the size of the VCF considerably.  The final VCF contains all SNP, indel or zero mapped/coverage positions
awk ' $0 ~ /#/ || $8 !~ /^AN=2;/ {print $0}' ${n}.allsites.vcf > ${n}.ready-mem.vcf
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${igvtools} index ${n}.ready-mem.vcf

echo "***Deleting Files"
rm -rf spades_output/
rm ${n}.unsortSNPsZeroCoverage.vcf
rm ${n}.all.bam
rm ${n}.dup.bam
rm ${n}.dup.bam.bai
rm ${n}.sorted.bam
rm ${n}.sorted.bam.bai
rm ${n}.realignedBam.bam
rm ${n}.realignedBam.bai
rm $forReads
rm $revReads
rm igv.log
rm ${n}.allsites.vcf
rm ${n}.allsites.vcf.idx
rm ${n}.forIndelRealigner.intervals
rm ${n}.recal_data.grp

rm ${n}.SNPsMapzeroNoHeader.vcf
rm ${n}.vcfFormated
rm ${n}.keepTheseZeroCovPositions
rm ${n}.duplicates
rm ${n}.zeroCoveragePositions
rm ${n}.calledSNPpositions
rm ${n}.body
rm ${n}.header
rm ${n}.hapreadyOnlySNPs.vcf

###################################
# The next 5 steps collect metrics
###################################

#Quality Score Distribution
echo "***Quality Score Distribution"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} QualityScoreDistribution REFERENCE_SEQUENCE=$ref INPUT=${n}.ready-mem.bam CHART_OUTPUT=${n}.QualityScorceDistribution.pdf OUTPUT=${n}.QualityScoreDistribution ASSUME_SORTED=true

#Mean Quality by Cycle
echo "***Mean Quality by Cycle"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} CollectMultipleMetrics REFERENCE_SEQUENCE=$ref INPUT=${n}.ready-mem.bam OUTPUT=${n}.Quality_by_cycle PROGRAM=MeanQualityByCycle ASSUME_SORTED=true

#Collect Alignment Summary Metrics
echo "***Collect Alignment Summary Metrics"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=$ref INPUT=${n}.ready-mem.bam OUTPUT=${n}.AlignmentMetrics ASSUME_SORTED=true

#Collect GC Bias Error
echo "***Collect GC Bias Error"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} CollectGcBiasMetrics REFERENCE_SEQUENCE=$ref INPUT=${n}.ready-mem.bam OUTPUT=${n}.CollectGcBiasMetrics CHART_OUTPUT=${n}.GC.PDF SUMMARY_OUTPUT=${n}.GC.summary.txt ASSUME_SORTED=true

#Collect Insert Size Metrics
echo "***Collect Insert Size Metrics"
$HOME/prog/jdk1.8.0_72/bin/./java -Xmx4g -jar ${picard} CollectInsertSizeMetrics REFERENCE_SEQUENCE=$ref INPUT=${n}.ready-mem.bam HISTOGRAM_FILE=${n}.InsertSize.pdf OUTPUT=${n}.CollectInsertSizeMetrics ASSUME_SORTED=true

cat ${n}.AlignmentMetrics > ${n}.Metrics_summary.txt
cat ${n}.CollectInsertSizeMetrics >> ${n}.Metrics_summary.txt

echo "***Organizing files"

#Move to qualityvalues subfolder
mkdir qualityvalues
mv ${n}.GC.PDF qualityvalues/
mv ${n}.QualityScorceDistribution.pdf qualityvalues/
mv ${n}.InsertSize.pdf qualityvalues/
mv ${n}.Quality_by_cycle*.pdf qualityvalues/

#Remove files
rm ${n}.CollectInsertSizeMetrics
rm ${n}.Quality_by_cycle.quality_distribution_metrics
rm ${n}.Quality_by_cycle.quality_by_cycle_metrics
rm ${n}.Quality_by_cycle.alignment_summary_metrics
rm ${n}.CollectGcBiasMetrics
rm ${n}.QualityScoreDistribution
rm ${n}.coverage

###########################
echo "***Getting stats for $n"

echo "fastq.gz file sizes:" > ${n}.stats2.txt
du -sh ../zips/*.gz | sed 's%../zips%%' >> ${n}.stats2.txt

echo "Unmapped fastq file sizes:" >> ${n}.stats2.txt
du -sh ../unmappedReads/*.gz  | sed 's%../unmappedReads%%' >> ${n}.stats2.txt

echo "Unmapped contig count:" >> ${n}.stats2.txt
grep -c ">" ./spades_output/scaffolds.fasta >> ${n}.stats2.txt
echo "" >> ${n}.stats2.txt

#some stats
sed -n 7,8p ${n}.Metrics_summary.txt | awk '{print $2}' >> ${n}.stats2.txt
sed -n 7,8p ${n}.Metrics_summary.txt | awk '{print $3}' >> ${n}.stats2.txt
sed -n 7,8p ${n}.Metrics_summary.txt | awk '{print $8}' >> ${n}.stats2.txt
readcount=$(sed -n 8p ${n}.Metrics_summary.txt | awk '{print $3}')
echo "" >> ${n}.stats2.txt

echo "***Bamtools running"
aveCoverage=$($HOME/prog/bamtools/bin/bamtools-2.4.0 coverage -in ${n}.ready-mem.bam | awk '{sum+=$3} END { print sum/NR"X"}')
echo "Average depth of coverage: $aveCoverage" >> ${n}.stats2.txt

#genome coverage
percGenomeMissing=$(awk -v x="$zeroposition" -v y="$refsize" 'BEGIN { print(x/y)*100}')
percGenomeCoverage="$(echo "100 - $percGenomeMissing" | bc)"
echo "Percent of reference with coverage: ${percGenomeCoverage}%" >> ${n}.stats2.txt

#cat ${n}.stats2.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" >> ${n}.stats.txt
cat ${n}.stats2.txt > ${n}.stats.txt
echo "" >> ${n}.stats.txt

rm ${n}.stats2.txt
###########################

#  Add Insert_Size and Read_Length to stats.txt file
echo 'Mean_Insert_Size  Standard_Deviation:' >> ${n}.stats.txt
awk 'BEGIN {OFS="\t"} { print $5,$6 }' ${n}.Quality_by_cycle.insert_size_metrics | awk 'FNR == 8 {print $0}' >> ${n}.stats.txt

echo 'Mean_Read_Length:' >> ${n}.stats.txt
awk 'BEGIN {OFS="\t"} { print $16 }' ${n}.AlignmentMetrics | awk 'FNR == 10 {print $0}' >> ${n}.stats.txt

echo "" >> ${n}.stats.txt

#  Add SNP call numbers to stats.txt file
echo "SNP and zero coverage positions in ${n}.SNPsZeroCoverage.vcf:" >> ${n}.stats.txt
egrep -v "#" ${n}.SNPsZeroCoverage.vcf | grep -c ".*" >> ${n}.stats.txt

echo "SNPs of AC2 and QUAL > 300:" >> ${n}.stats.txt
egrep -v "#" ${n}.SNPsZeroCoverage.vcf | egrep "AC=2" | awk '$6 > 300' | grep -c ".*" >> ${n}.stats.txt

#  Show Mean Coverage at Terminal and coverageReport
echo "Mean Coverage"

echo "Sample identified and ran as:  $1" >> $HOME/pipeline/public-master/report/dailyReport.txt

echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% "
echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% " >> $HOME/pipeline/public-master/report/coverageReport.txt
echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% " >> $HOME/pipeline/public-master/report/dailyReport.txt

mv ${n}.Metrics_summary.txt qualityvalues/
mv ${n}.stats.txt qualityvalues/
rm ${n}.Quality_by_cycle.insert_size_metrics
rm ${n}.AlignmentMetrics

if [ -e ../*identifier_out.txt ];then
    cat ../*out1.txt ../*out2.txt > ../${n}-identification.txt
    rm ../*identifier_out1*
    rm ../*identifier_out2*
fi

mv ${startingdir}/fastq ${startingdir}/spoligo
rm ${startingdir}/spoligo/*fastq
rm -r ${startingdir}/temp
ln qualityvalues/${n}.stats.txt ./stats-${n}.txt

cp $0 ./

#Make dailyStats.txt for each stats.txt made for each isolate.
echo "" >> $HOME/pipeline/public-master/report/dailyStats.txt
echo "" >> $HOME/pipeline/public-master/report/dailyStats.txt
echo "" >> $HOME/pipeline/public-master/report/dailyStats.txt
echo "ADD_MARKER" >> $HOME/pipeline/public-master/report/dailyStats.txt
echo "" >> $HOME/pipeline/public-master/report/dailyStats.txt
echo "<------- $n $1 ------->" >> $HOME/pipeline/public-master/report/dailyStats.txt
cat qualityvalues/${n}.stats.txt >> $HOME/pipeline/public-master/report/dailyStats.txt

echo "**************************** END $n ****************************"

#
#  Created by Stuber, Tod P - APHIS on 11/08/12.
#
