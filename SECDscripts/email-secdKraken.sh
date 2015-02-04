#! /bin/sh


echo "Start time: `date`" > /scratch/report/dailyTime 
starttime=`date +%s`

##############################
# Environment controls:

if [[ $1 == all ]]; then  #All Viruses
	organism="all"

elif [[ $1 == secd ]]; then #4 porcine coronaviruses
	organism="secd"

else
	echo ""
	echo "Incorrect argument! Must use one of the following arguments: all, secd"
	echo: "For example, type ~$ email-kaityKraken.sh secd"
	echo ""
	exit 1
fi

if [[ $2 == "brien" ]]; then #Email to Kaity only
    email_list="kaitlin.e.brien@aphis.usda.gov"
elif [[ $2 == "stuber" ]]; then #Email to Tod only
    email_list="tod.p.stuber@aphis.usda.gov"
elif [[ $2 == "stuber+brien" ]]; then
    email_list="kaitlin.e.brien@aphis.usda.gov, tod.p.stuber@aphis.usda.gov"
elif [[ $2 == "all" ]]; then 
    email_list="kaitlin.e.brien@aphis.usda.gov, tod.p.stuber@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"
else
    echo ""
    echo "Incorrect argument! Must provide an argument to email report: brien, stuber, stuber+brien, all"
    echo "For example: email-secdKraken.sh secd all"
    echo ""
    exit 1
fi


function loopfiles () {

currentdir=`pwd`

count='0'
for i in *.fastq*; do
	n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
	echo "n is: $n"
	mkdir -p $n
	mv $i $n/
done

for f in *; do
	cd $currentdir
	cd $f
	echo "In directory $f"
	constrainedKraken.sh $organism &
done

}
rm /scratch/report/kraken/* 
loopfiles &&
wait

echo "End Time: `date`" >> /scratch/report/dailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "Run time: $runtime seconds" >> /scratch/report/dailyTime
echo "" >> /scratch/report/dailyTime

cd $currentdir
touch sampleSummary.txt
echo "******************** Summary of Samples ********************" >> sampleSummary.txt
echo "" >> sampleSummary.txt
for i in `find . -name "reportHeader.txt"`; do
    cat $i >> sampleSummary.txt
    rm $i
done  

echo "******************** Detailed Report ********************" >> sampleSummary.txt
echo "" >> sampleSummary.txt
for i in `find . -name "*resultsSummary.txt"`; do
    cat $i >> sampleSummary.txt
done

find . -name "*.pdf" > list
find . -name "*KronaGraphic.html" >> list

cat /scratch/report/dailyTime > /scratch/report/kraken/email-report.txt
cat $currentdir/sampleSummary.txt >> /scratch/report/kraken/email-report.txt

#email_list="kaitlin.e.brien@aphis.usda.gov tod.p.stuber@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"
#email_list="kaitlin.e.brien@aphis.usda.gov"
#email_list="kaitlin.e.brien@aphis.usda.gov tod.p.stuber@aphis.usda.gov"

cat /scratch/report/kraken/email-report.txt | mutt -s "Sequence submission results" -a `cat list` -- $email_list

rm list


# Created by Kaity Brien, 2014-11-14
