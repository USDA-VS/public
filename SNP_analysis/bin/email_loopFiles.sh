#!/bin/sh

#  email_loopFiles.sh

echo "Start Time: `date`" > /scratch/report/dailyTime
starttime=`date +%s`

`/home/tstuber/workspace/public/SNP_analysis/bin/loopFiles.sh` &&

echo "" >> /scratch/report/dailyTime
echo "End Time: `date`" >> /scratch/report/dailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "Run time: $runtime seconds" >> /scratch/report/dailyTime

echo "e-mailing files"

cat /scratch/report/dailyTime > /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt
echo "ADD_MARKER" >> /scratch/report/email_processZips.txt
echo "" >> /scratch/report/dailyReport.txt
cat /scratch/report/dailyReport.txt >> /scratch/report/email_processZips.txt

cat /scratch/report/spoligoCheck.txt >> /scratch/report/email_processZips.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/email_processZips.txt

echo "ADD_MARKER" >> /scratch/report/email_processZips.txt

cat /scratch/report/dailyStats.txt >> /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt

grep -v '*' /scratch/report/email_processZips.txt | grep -v "Stats for BAM file" | sed 's/ADD_MARKER/******************************************/g' > /scratch/report/email_processZips2.txt

email_list="tod.p.stuber@aphis.usda.gov"
#email_list="tod.p.stuber@aphis.usda.gov patrick.m.camp@aphis.usda.gov David.T.Farrell@aphis.usda.gov Christine.R.Quance@aphis.usda.gov suelee.robbe-austerman@aphis.usda.gov"

cat /scratch/report/email_processZips2.txt | mutt -s "WGS results" -- $email_list

#mail -s "WGS results" tod.p.stuber@aphis.usda.gov < /scratch/report/email_processZips2.txt

date >> /scratch/report/mlstCheck_all.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/mlstCheck_all.txt

#echo "<html>" > /scratch/report/assembly_stats.html
#awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' /scratch/report/email_processZips2.txt >> /scratch/report/assembly_stats.html
#echo "</html>" >> /scratch/report/assembly_stats.html

# E-mail results
#basen=`basename $0`
# As attachment
#echo "See attachment"| mutt -s "WGS results" -a /scratch/report/assembly_stats.html -- tod.p.stuber@aphis.usda.gov

#echo "See attachment"| mutt -s "WGS results" -a /scratch/report/assembly_stats.html -- patrick.m.camp@aphis.usda.gov

#echo "See attachment"| mutt -s "WGS results" -a /scratch/report/assembly_stats.html -- suelee.robbe-austerman@aphis.usda.gov

#echo "See attachment"| mutt -s "WGS results" -a /scratch/report/assembly_stats.html -- Christine.R.Quance@aphis.usda.gov

#echo "See attachment"| mutt -s "WGS results" -a /scratch/report/assembly_stats.html -- David.T.Farrell@aphis.usda.gov


#
#  Created by Tod Stuber on 11/09/12.
#
