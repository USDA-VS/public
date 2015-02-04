#!/usr/bin/env python

import xml.etree.ElementTree as ET
from sys import argv

script, input = argv

mytree = ET.parse(input)
myroot = mytree.getroot()

for Blast_iteration in myroot.findall('BlastOutput_iterations'):
    for Iteration in Blast_iteration.findall('Iteration'):
        queryID = Iteration.find('Iteration_query-def').text.split("_")
        #length = queryID[-3]
        #coverage = queryID[-1]
        #totalseq = int(length) * float(coverage)
        for Iteration_hits in Iteration.findall('Iteration_hits'):
            for Hit in Iteration_hits.findall('Hit'):
                Hit_def = Hit.find('Hit_def').text
                Hit_def_remove = Hit_def.replace("PREDICTED: ", "")
                Hit_def_split = Hit_def_remove.split(' ')
                Hit_accession = Hit.find('Hit_accession').text
                Hit_len = Hit.find('Hit_len').text
                for Hit_hsps in Hit.findall('Hit_hsps'):
                    for Hsp in Hit_hsps.findall('Hsp'):
                        Hsp_bitscore = Hsp.find('Hsp_bit-score').text
                        Hsp_evalue = Hsp.find('Hsp_evalue').text
                print queryID, Hit_def_split[0], Hit_def_split[1], Hit_accession, Hit_len, Hsp_bitscore, Hsp_evalue, Hit_def
