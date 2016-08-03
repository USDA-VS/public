#!/bin/sh

for f in *.fasta; do 
grep -v "^$" $f > file
h=`grep -c ">" file`
a=`grep -c ".*" file`
echo "h: $h a: $a"
singlelinetest=`expr $a / $h`
echo "singlelingtest: $singlelinetest"

if [[ $singlelinetest == 2 ]]; then
        echo "contig file ( $f ) is formatted correctly"
        cat $f > ${n}-contigsoriginal.fa
else
        # If contigs going in are not all on the same line then newlines are removed to put fasta on a single line, but the name header name is changed.
        echo "contig file ( $f ) was reformated placing sequence on a sigle line"
        awk '{if ($0 ~ />/ ) {print ">"} else print $0}' $1 | tr -d "\n" | sed -e 's:>:\n>\n:g' | awk '{if ($0 ~  />/ ) {print ">contigs-" x++} else print $0}' | grep -v "^$" > ${n}-contigsoriginal.fa
fi
rm file


# created 2015-08-17
