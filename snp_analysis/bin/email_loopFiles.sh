#!/bin/sh

#  email_loopFiles.sh

echo "Start Time: $(date)" > $HOME/pipeline/public-master/report/dailyTime
starttime=$(date +%s)

echo "Please wait.  Searching for TB complex, Brucella and paratuberculosis oligos and then starting appropriate processZips.sh argument"
`loopFiles.sh` &&

echo "" >> $HOME/pipeline/public-master/report/dailyTime
echo "End Time: $(date)" >> $HOME/pipeline/public-master/report/dailyTime
endtime=$(date +%s)
runtime=$(($endtime - $starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> $HOME/pipeline/public-master/report/dailyTime

echo "e-mailing files"

cat $HOME/pipeline/public-master/report/dailyTime > $HOME/pipeline/public-master/report/email_processZips.txt
echo "" >> $HOME/pipeline/public-master/report/email_processZips.txt
echo "ADD_MARKER" >> $HOME/pipeline/public-master/report/email_processZips.txt
echo "" >> $HOME/pipeline/public-master/report/dailyReport.txt
cat $HOME/pipeline/public-master/report/dailyReport.txt >> $HOME/pipeline/public-master/report/email_processZips.txt

cat $HOME/pipeline/public-master/report/spoligoCheck.txt >> $HOME/pipeline/public-master/report/email_processZips.txt
cat $HOME/pipeline/public-master/report/mlstCheck.txt >> $HOME/pipeline/public-master/report/email_processZips.txt

echo "ADD_MARKER" >> $HOME/pipeline/public-master/report/email_processZips.txt

cat $HOME/pipeline/public-master/report/dailyStats.txt >> $HOME/pipeline/public-master/report/email_processZips.txt
echo "" >> $HOME/pipeline/public-master/report/email_processZips.txt

grep -v "*" $HOME/pipeline/public-master/report/email_processZips.txt | grep -v "Stats for BAM file" \
| sed 's/ADD_MARKER/******************************************/g' > $HOME/pipeline/public-master/report/email_processZips2.txt

if [[ $1 == me ]]; then
	email_list="marc-olivier.duceppe@inspection.gc.ca"
	else
	email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaea@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"
fi

cat $HOME/pipeline/public-master/report/email_processZips2.txt | mail -s "WGS results" $email_list

date >> $HOME/pipeline/public-master/report/mlstCheck_all.txt
cat $HOME/pipeline/public-master/report/mlstCheck.txt >> $HOME/pipeline/public-master/report/mlstCheck_all.txt


#
#  Created by Tod Stuber on 11/09/12.
#
