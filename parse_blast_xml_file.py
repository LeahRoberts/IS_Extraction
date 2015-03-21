#!/usr/bin/env python

############## Python script to parse the xml output of a blastn comparison #####################
# Author: Leah Roberts
# Affiliation: Beatson lab group, University of Queensland St Lucia
# Date: March 2015

# Script Description:
# Simple script for parsing out the name of a blastn comparison along with the start and end coordinates
# of the blastn match.

# Usage:
# To generate the blastn xml files, see blast commandline documentation.
# To run the script on one xml file, use the command:
# $ python parse_blast_xml_file.py <file_name.xml>

# If running the script on multiple xml files, simply run a loop:
# $ for f in *; do python parse_blast_xml_file.py $f; done

# The output of the script will be text files containing the query, start and end coordinates,  and the
# evalue, eg:

# Query:IS1.fna.results.xml
# sequence:gnl|BL_ORD_ID|0 EMBOSS_001.1
# 6013..6780
# Evalue:0.0


from Bio.Blast import NCBIXML
import sys

name = sys.argv[1].split(".")[0]
blast_xml = sys.argv[1]

file = open(name + ".txt", "w")
blast_records = NCBIXML.parse(open(blast_xml, "r"))
for blast_record in blast_records:
        for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                        file.write("Query:" + sys.argv[1] + "\n" + "sequence:" + alignment.title + "\n" + "Start:" + str(hsp.sbjct_start) + "\n" + "End:" + str(hsp.sbjct_end) + "\n" + "Evalue:" + str(hsp.expect) + "\n")
file.close()

exit()
