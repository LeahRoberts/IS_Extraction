#!/usr/bin/env python

# Author: Leah Roberts
# Affiliation: Scott Beatson Lab Group University of Queensland St Lucia
# Date: April 2014

####################### Remove IS - keep flanking region ########################

###### This script is designed to parse results from "Extract IS" script ########

### Description:

# This script takes the output of the "Extract IS" script, i.e. a fasta file with
# an IS sequence and 1000 bp flanking region either side, and parses out just the
# flanking region without the IS sequence.

### Usage:

# Run the "Extract_IS" script first, which will generate a 'results' folder with the
# IS and the flanking regions. 

# Run this script in the directory above the 'results' folder:
# $ python Remove_IS_keep_flank.py

# The output isn't written to a file, so will be printed to the screen.
# To write to a file, parse the output at the end of the command:
# $ python Remove_IS_keep_flank.py > flanking_regions.fa


import Bio
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Takes each of the IS files from a 'results' directory. Each IS is in a separate .fna file
# Can print the in_file name to test that the script is working

flank = 100
for in_file in glob.glob('results/*'):
#	print in_file

# Takes the string and parses the header and the fasta sequence.
# fixes the string twice to firstly start from the 1000th position, and then end 1000 positions from the end.
# Prints the results for each strain which can be parsed into a concatenated file.

	for record in SeqIO.parse(in_file, "fasta"):
		cur = str(record.seq)
		ID = record.id
		leftflank = cur[0:flank]
		rightflank = cur[-flank:]
		total = leftflank + rightflank

		print ">" + ID
		print total
		
# The header may need to be fixed with a simple sed command (the headers had the latter part missing)
