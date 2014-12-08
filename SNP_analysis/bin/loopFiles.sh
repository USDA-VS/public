#/bin/sh

#  bovis_loopFiles.sh
#  Working directory must have Brucella or Bovis zip files
#  Brucelle files must begin with a "B"
#  Bovis files must NOT begin with a "B"


date >> /scratch/report/coverageReport.txt
echo "" > /scratch/report/dailyReport.txt

# Reset spoligo and bruc mlst check file
echo "" > /scratch/report/spoligoCheck.txt
echo "" > /scratch/report/mlstCheck.txt
echo "WG Spoligo Check" >> /scratch/report/spoligoCheck.txt
echo "Brucella MLST Check" >> /scratch/report/mlstCheck.txt

#Reset file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "TB Number" "Octal Code" > "/bioinfo11/TStuber/Results/_Mycobacterium/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

echo "Isolate Total_Bases AveDep %>Q15" | awk '{ printf("%-12s %-12s %-10s %-10s\n", $1, $2, $3, $4) }' >> /scratch/report/dailyReport.txt
echo "" >> /scratch/report/dailyReport.txt
echo ""  > /scratch/report/dailyStats.txt

currentdir=`pwd`

for i in *.fastq.gz; do

n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
echo "n is : $n"

mkdir -p $n
mv $i $n/
done

for f in *; do
echo "The cd is $currentdir"
cd $currentdir
t=`echo $f | sed 's/\(^.\).*/\1/'`
echo $f


if [ $t == B ]

#  Files beginning with a B are process as Brucella
	then

	echo "***$f has been started as Brucella"
        
        
###
    cd ./$f
    mkdir ./temp
    cp *R1*.fastq.gz ./temp

    `gunzip ./temp/*R1*.fastq.gz && bruc_oligo_identifier.sh ./temp/*R1*.fastq | tee tee_bruc_oligo_identifier_out1.txt` &

#`cd ./$f; /Users/Shared/_programs/_my_scripts/bruc_oligo_identifier.sh $R1; cd ..` &

###

#`cd ./$f; /Users/Shared/_programs/_my_scripts/brucsuis_processZips.sh; cd ..` &
#elif [ $t == T ]

#  Files beginning with a T are process as Taylorella
#	then

#	echo "***$f has been started as Taylorella"

#`cd ./$f; processZips.sh taylorella; cd ..` &

elif [ $t == P ]

#  Files beginning with a P are process as Pastuerella
then

echo "***$f has been started as Pastuerella"

`cd ./$f; /Users/Shared/_programs/_my_scripts/past_processZips.sh; cd ..` &


#elif [ $t == C ]

#  Files beginning with a P are process as Pastuerella
#then

#echo "***$f has been started as Pastuerella"

#`cd ./$f; /Users/Shared/_programs/_my_scripts/mem_campy_processZips.sh; cd ..` &

#  Files NOT beginning with a B are process as Bovis

else

echo "***$f has been started as Bovis"


`cd ./$f; processZips.sh bovis; cd ..` &

#`cd ./$f; /Users/Shared/_programs/_my_scripts/bovis_processZips.sh; cd..` &

fi

done


#
#  Created by Tod Stuber on 11/05/12.
#
