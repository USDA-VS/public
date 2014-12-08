#!/bin/sh

#################################################################################
#  Dependencies ---
#   Unix..ish setup :)
#################################################################################

echo "**********************START**********************"

spacer01="TGATCCAGAGCCGGCGACCCTCTAT|ATAGAGGGTCGCCGGCTCTGGATCA"
spacer02="CAAAAGCTGTCGCCCAA|TTGGGCGACAGCTTTTG"
spacer03="CCGTGCTTCCAGTGATCGCCTTCTA|TAGAAGGCGATCACTGGAAGCACGG"
spacer04="ACGTCATACGCCGACCAATCATCAG|CTGATGATTGGTCGGCGTATGACGT"
spacer05="TTTTCTGACCACTTGTGCGGGATTA|TAATCCCGCACAAGTGGTCAGAAAA"
spacer06="CGTCGTCATTTCCGGCTTCAATTTC|GAAATTGAAGCCGGAAATGACGACG"
spacer07="GAGGAGAGCGAGTACTCGGGGCTGC|GCAGCCCCGAGTACTCGCTCTCCTC"
spacer08="CGTGAAACCGCCCCCAGCCTCGCCG|CGGCGAGGCTGGGGGCGGTTTCACG"
spacer09="ACTCGGAATCCCATGTGCTGACAGC|GCTGTCAGCACATGGGATTCCGAGT"
spacer10="TCGACACCCGCTCTAGTTGACTTCC|GGAAGTCAACTAGAGCGGGTGTCGA"
spacer11="GTGAGCAACGGCGGCGGCAACCTGG|CCAGGTTGCCGCCGCCGTTGCTCAC"
spacer12="ATATCTGCTGCCCGCCCGGGGAGAT|ATCTCCCCGGGCGGGCAGCAGATAT"
spacer13="GACCATCATTGCCATTCCCTCTCCC|GGGAGAGGGAATGGCAATGATGGTC"
spacer14="GGTGTGATGCGGATGGTCGGCTCGG|CCGAGCCGACCATCCGCATCACACC"
spacer15="CTTGAATAACGCGCAGTGAATTTCG|CGAAATTCACTGCGCGTTATTCAAG"
spacer16="CGAGTTCCCGTCAGCGTCGTAAATC|GATTTACGACGCTGACGGGAACTCG"
spacer17="GCGCCGGCCCGCGCGGATGACTCCG|CGGAGTCATCCGCGCGGGCCGGCGC"
spacer18="CATGGACCCGGGCGAGCTGCAGATG|CATCTGCAGCTCGCCCGGGTCCATG"
spacer19="TAACTGGCTTGGCGCTGATCCTGGT|ACCAGGATCAGCGCCAAGCCAGTTA"
spacer20="TTGACCTCGCCAGGAGAGAAGATCA|TGATCTTCTCTCCTGGCGAGGTCAA"
spacer21="TCGATGTCGATGTCCCAATCGTCGA|TCGACGATTGGGACATCGACATCGA"
spacer22="ACCGCAGACGGCACGATTGAGACAA|TTGTCTCAATCGTGCCGTCTGCGGT"
spacer23="AGCATCGCTGATGCGGTCCAGCTCG|CGAGCTGGACCGCATCAGCGATGCT"
spacer24="CCGCCTGCTGGGTGAGACGTGCTCG|CGAGCACGTCTCACCCAGCAGGCGG"
spacer25="GATCAGCGACCACCGCACCCTGTCA|TGACAGGGTGCGGTGGTCGCTGATC"
spacer26="CTTCAGCACCACCATCATCCGGCGC|GCGCCGGATGATGGTGGTGCTGAAG"
spacer27="GGATTCGTGATCTCTTCCCGCGGAT|ATCCGCGGGAAGAGATCACGAATCC"
spacer28="TGCCCCGGCGTTTAGCGATCACAAC|GTTGTGATCGCTAAACGCCGGGGCA"
spacer29="AAATACAGGCTCCACGACACGACCA|TGGTCGTGTCGTGGAGCCTGTATTT"
spacer30="GGTTGCCCCGCGCCCTTTTCCAGCC|GGCTGGAAAAGGGCGCGGGGCAACC"
spacer31="TCAGACAGGTTCGCGTCGATCAAGT|ACTTGATCGACGCGAACCTGTCTGA"
spacer32="GACCAAATAGGTATCGGCGTGTTCA|TGAACACGCCGATACCTATTTGGTC"
spacer33="GACATGACGGCGGTGCCGCACTTGA|TCAAGTGCGGCACCGCCGTCATGTC"
spacer34="AAGTCACCTCGCCCACACCGTCGAA|TTCGACGGTGTGGGCGAGGTGACTT"
spacer35="TCCGTACGCTCGAAACGCTTCCAAC|GTTGGAAGCGTTTCGAGCGTACGGA"
spacer36="CGAAATCCAGCACCACATCCGCAGC|GCTGCGGATGTGGTGCTGGATTTCG"
spacer37="CGCGAACTCGTCCACAGTCCCCCTT|AAGGGGGACTGTGGACGAGTTCGCG"
spacer38="CGTGGATGGCGGATGCGTTGTGCGC|GCGCACAACGCATCCGCCATCCACG"
spacer39="GACGATGGCCAGTAAATCGGCGTGG|CCACGCCGATTTACTGGCCATCGTC"
spacer40="CGCCATCTGTGCCTCATACAGGTCC|GGACCTGTATGAGGCACAGATGGCG"
spacer41="GGAGCTTTCCGGCTTCTATCAGGTA|TACCTGATAGAAGCCGGAAAGCTCC"
spacer42="ATGGTGGGACATGGACGAGCGCGAC|GTCGCGCTCGTCCATGTCCCACCAT"
spacer43="CGCAGAATCGCACCGGGTGCGGGAG|CTCCCGCACCCGGTGCGATTCTGCG"

# Starting working directory must be BWA-GATK folder with included /Zips file containing 2 zipped fastq files.
echo "directory"
pwd
# Make fastq directory
mkdir ./../fastq
cp ./../Zips/*.fastq.gz ./../fastq

# change working directory to /fastq
cd ./../fastq

echo "starting to unzip files"
# Unzip files
gunzip *.fastq.gz
echo "finished unzipping files"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'` #grab name minus the .vcf

sp01=`egrep $spacer01 $forReads $revReads | wc -l`
sp02=`egrep $spacer02 $forReads $revReads | wc -l`
sp03=`egrep $spacer03 $forReads $revReads | wc -l`
sp04=`egrep $spacer04 $forReads $revReads | wc -l`
sp05=`egrep $spacer05 $forReads $revReads | wc -l`
sp06=`egrep $spacer06 $forReads $revReads | wc -l`
sp07=`egrep $spacer07 $forReads $revReads | wc -l`
sp08=`egrep $spacer08 $forReads $revReads | wc -l`
sp09=`egrep $spacer09 $forReads $revReads | wc -l`
sp10=`egrep $spacer10 $forReads $revReads | wc -l`
sp11=`egrep $spacer11 $forReads $revReads | wc -l`
sp12=`egrep $spacer12 $forReads $revReads | wc -l`
sp13=`egrep $spacer13 $forReads $revReads | wc -l`
sp14=`egrep $spacer14 $forReads $revReads | wc -l`
sp15=`egrep $spacer15 $forReads $revReads | wc -l`
sp16=`egrep $spacer16 $forReads $revReads | wc -l`
sp17=`egrep $spacer17 $forReads $revReads | wc -l`
sp18=`egrep $spacer18 $forReads $revReads | wc -l`
sp19=`egrep $spacer19 $forReads $revReads | wc -l`
sp20=`egrep $spacer20 $forReads $revReads | wc -l`
sp21=`egrep $spacer21 $forReads $revReads | wc -l`
sp22=`egrep $spacer22 $forReads $revReads | wc -l`
sp23=`egrep $spacer23 $forReads $revReads | wc -l`
sp24=`egrep $spacer24 $forReads $revReads | wc -l`
sp25=`egrep $spacer25 $forReads $revReads | wc -l`
sp26=`egrep $spacer26 $forReads $revReads | wc -l`
sp27=`egrep $spacer27 $forReads $revReads | wc -l`
sp28=`egrep $spacer28 $forReads $revReads | wc -l`
sp29=`egrep $spacer29 $forReads $revReads | wc -l`
sp30=`egrep $spacer30 $forReads $revReads | wc -l`
sp31=`egrep $spacer31 $forReads $revReads | wc -l`
sp32=`egrep $spacer32 $forReads $revReads | wc -l`
sp33=`egrep $spacer33 $forReads $revReads | wc -l`
sp34=`egrep $spacer34 $forReads $revReads | wc -l`
sp35=`egrep $spacer35 $forReads $revReads | wc -l`
sp36=`egrep $spacer36 $forReads $revReads | wc -l`
sp37=`egrep $spacer37 $forReads $revReads | wc -l`
sp38=`egrep $spacer38 $forReads $revReads | wc -l`
sp39=`egrep $spacer39 $forReads $revReads | wc -l`
sp40=`egrep $spacer40 $forReads $revReads | wc -l`
sp41=`egrep $spacer41 $forReads $revReads | wc -l`
sp42=`egrep $spacer42 $forReads $revReads | wc -l`
sp43=`egrep $spacer43 $forReads $revReads | wc -l`

echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43" > $n.spacer.txt

echo "$sp01	$sp02	$sp03	$sp04	$sp05	$sp06	$sp07	$sp08	$sp09	$sp10	$sp11	$sp12	$sp13	$sp14	$sp15	$sp16	$sp17	$sp18	$sp19	$sp20	$sp21	$sp22	$sp23	$sp24	$sp25	$sp26	$sp27	$sp28	$sp29	$sp30	$sp31	$sp32	$sp33	$sp34	$sp35	$sp36	$sp37	$sp38	$sp39	$sp40	$sp41	$sp42	$sp43" >> $n.spacer.txt

cat $n.spacer.txt | awk 'NR==2 {for(i=1;i<=NF;i++) if ($i >= 5) print 1; else print 0}' | tr -cd "[:print:]" | fold -w3 > $n.myspacers

chmod 755 ./$n.myspacers

mybinaries=`cat ./$n.myspacers`

for i in $mybinaries; do 
if [ $i == 000 ]
then
	echo "0" >> $n.octalcode.txt
elif [ $i == 001 ]
then
	echo "1" >> $n.octalcode.txt
elif [ $i == 010 ]
	then
	echo "2" >> $n.octalcode.txt
elif [ $i == 011 ]
	then
	echo "3" >> $n.octalcode.txt
elif [ $i == 100 ]
	then
	echo "4" >> $n.octalcode.txt
elif [ $i == 101 ]
	then
	echo "5" >> $n.octalcode.txt
elif [ $i == 110 ]
	then
	echo "6" >> $n.octalcode.txt	
elif [ $i == 111 ]
	then	
	echo "7" >> $n.octalcode.txt
elif [ $i == 0 ]
	then	
	echo "0" >> $n.octalcode.txt
elif [ $i == 1 ]
	then
# Changed 4 to 1 on 2013-10-25
	echo "1" >> $n.octalcode.txt
else
	echo "***Error***" >> $n.octalcode.txt
fi
done

WGSpoligo=`cat $n.octalcode.txt | tr -cd "[:print:]"`

# Grab blot spoligo Octal code
#p=`grep "$n" /Volumes/Data_HD/Mycobacterium/bovis_tags.txt`
#fileName=`echo $p | awk '{print $4}'`
#fileName2=`echo $fileName | tr -d "[:space:]"`

# Add infor to spoligoCheck.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck.txt
#echo "Blot:		$fileName2" >> /scratch/report/spoligoCheck.txt

#Make FileMaker file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "$n" "$WGSpoligo" >> "/bioinfo11/TStuber/Results/_Mycobacterium/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

#if [ $WGSpoligo == $fileName2 ]
#then
#	echo "PASSED" >> /scratch/report/spoligoCheck.txt
#else
#	echo "FAILED" >> /scratch/report/spoligoCheck.txt
#fi
#echo "" >> /scratch/report/spoligoCheck.txt

# Add infor to spoligoCheck_all.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck_all.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck_all.txt
#echo "Blot:		$fileName2" >> /scratch/report/spoligoCheck_all.txt

#if [ $WGSpoligo == $fileName2 ]
#then
#	echo "PASSED" >> /scratch/report/spoligoCheck_all.txt
#else
#	echo "FAILED" >> /scratch/report/spoligoCheck_all.txt
#fi
#echo "" >> /scratch/report/spoligoCheck_all.txt

# move back a directory to main sample folder
cd ..

#
#  Created by Stuber, Tod P - APHIS on 03/07/2013.
#
