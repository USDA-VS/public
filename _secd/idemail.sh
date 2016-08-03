#!/bin/sh

# Script will send out a summary email of stats output from idvirus.sh

# flag -m will email just "M"e
# flag -k will run Kraken
mflag=
kflag=
while getopts 'mk' OPTION; do
    case $OPTION in
        m) mflag=1
        ;;
        k) kflag=1
        ;;
        ?) echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))


if [ -z $1 ]; then 
    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: gen, testflu, allflu, sivall, h5n2, h5n8, h11n9, secd, reo, vsv, isav"
    echo ""
    echo "Set optional flags"
    echo -e '   flag -m will email just "M"e'
    echo -e '   flag -k wll run "K"raken'
    echo ""
    echo "Example: [prompt]$ idemail.sh aiall"
    echo ""	
    exit 1
fi

#####
argUsed=`echo $1 | tr '[:lower:]' '[:upper:]'`

# Clear email file
echo "Reference set: $argUsed" > /scratch/report/idemailsummary
echo "" >> /scratch/report/idemailsummary

echo "Start Time: `date`" > /scratch/report/iddailyTime
starttime=`date +%s`

for i in *.fastq*; do 
	n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
	echo "n is : $n"
	mkdir -p $n
	mv $i $n/
 done

if [ "$kflag" ]; then
	echo "Kraken ran"
	currentdir=`pwd`
	for f in *; do 
		cd $currentdir
		echo $f
		cd ./$f 
		idvirus.sh -ek $1 
	done
else
	echo "Kraken not ran"
        currentdir=`pwd`
        for f in *; do 
                cd $currentdir
                echo $f
                cd ./$f 
                idvirus.sh -e $1 
        done
fi


echo "End Time: `date`" >> /scratch/report/iddailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> /scratch/report/iddailyTime

cat /scratch/report/iddailyTime > /scratch/report/tempfile
echo "" >> /scratch/report/tempfile
cat /scratch/report/idemailsummary >> /scratch/report/tempfile
mv /scratch/report/tempfile /scratch/report/idemailsummary

if [ "$mflag" ]; then
	# just email me the summary
	email_list="tod.p.stuber@usda.gov"
	enscript /scratch/report/idemailsummary -B -j -r -f "Courier5" -o - | ps2pdf - /scratch/report/summary-report.pdf
        cat /scratch/report/idemailsummary | mutt -s "Samples analyzed" -a /scratch/report/summary-report.pdf -- $email_list
	
else
	# email all the summary
        email_list="tod.p.stuber@usda.gov Mary.L.Killian@aphis.usda.gov mia.kim.torchetti@aphis.usda.gov Suelee.Robbe-Austerman@aphis.usda.gov"
	enscript /scratch/report/idemailsummary -B -j -r -f "Courier5" -o - | ps2pdf - /scratch/report/summary-report.pdf
        cat /scratch/report/idemailsummary | mutt -s "Samples analyzed" -a /scratch/report/summary-report.pdf -- $email_list
fi


# created 2015-08-19 stuber
