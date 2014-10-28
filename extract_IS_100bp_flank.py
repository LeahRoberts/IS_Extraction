#!/usr/bin/env python

# Author: Mitchell Stanton-Cook (https://github.com/mscook)
# Affiliation: Scott Beatson Lab Group University of Queensland St Lucia
# Date: March 2014

############## Script to parse IS from embl file with flanking sequence ###################

### Description:

# This script takes in an embl file from a 'genomes' folder and parses out gene features,
# as well as flanking sequence of a predefined length, that match a specific feature type. 

## Usage:

# Execute script from above a 'genomes' directory containing embl files that you want to extract
# particular sequences from. 
# The header has been modified to suit usage with SeqFindR (https://github.com/mscook/seqfindr). 
# The amount of flanking sequence desired, as well as the feature type that you want extracted,
# can be changed to suit the user.


import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import os

in_files = glob.glob('genomes/*.embl')
flanking_region = 100
try:
        os.mkdir("results")
except OSError:
        print "a 'results' dir already exists"
        print "Overwriting"
stored = {}
for f in in_files:
        cur_genome = SeqIO.parse(f, "embl")
#print cur_genome
for record in cur_genome:
        for feat in record.features:
                if feat.type == 'mobile_element':
                        s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                        header = '>'+feat.qualifiers['mobile_element_type'][0].split(':')[-1]+","+feat.qualifiers['mobile_element_type'][0].split(':')[-1]+".."+str(s+1)+".."+str(e)+"("+str(strand)+"),""100bp flanked,[EC958 IS]"
                        flanked = FeatureLocation(s-flanking_region, e+flanking_region, strand)
                        out_seq = flanked.extract(record.seq)
                        fname = header[1:].split(',')[0].replace('unclassified','unc').replace('family', 'fam').replace('(', '').replace('partial', 'p').replace(')', '').replace(' ', '_').replace('/', '-').strip()+'.fna'
                        if fname in stored.keys():
                                old = fname
                                fname = fname.replace(".fna", "_"+str(stored[fname])+".fna")
                                stored[old] = stored[old]+1
                        else:
                                stored[fname] = 1
                        with open(os.path.join('results', fname), 'w') as out:
                                out.write(header+'\n')
                                out.write(str(out_seq)+'\n')
