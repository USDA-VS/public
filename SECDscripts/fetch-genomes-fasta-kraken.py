#!/usr/bin/python

import urllib2
import os
import sys
import time

if len(sys.argv) != 2:
    print "USAGE: fetch_genome.py <accession>"
    sys.exit(1)

url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text"

#os.mkdir(sys.argv[2])

id = sys.argv[1]
#for id in open(sys.argv[1]):
id = id.strip()
#if id == "":
#    continue

sys.stdout.write("Fetching %s..." % id)
sys.stdout.flush()
gbk_out_file = id + ".fasta"
#if os.path.exists(gbk_out_file):
#    print "already fetched"

open(gbk_out_file, "w").write(urllib2.urlopen(url_template % id).read())
print "Done"
time.sleep(1.0/3)

