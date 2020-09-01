#!/usr/bin/env python
# script ReGAl (Removal of Gaps in Alignments)written by Chandra Earl (sunray1@ufl.edu) to remove gaps in alignments based on a reference sequence

import os, sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from argparse import ArgumentParser

argp = ArgumentParser(description='Removes sites from an alignment where the reference sequence is empty')
argp.add_argument('-i', '--inalignment', help='name of the input alignment', required=True)
argp.add_argument('-f', '--fileformat', help='type of file input (ie. phylip-relaxed, fasta, clustal, stockholm)')
argp.set_defaults(fileformat='fasta')

args = argp.parse_args()
fileformat = args.fileformat
infile = args.inalignment



if fileformat == "phylip":
    fileformat = "phylip-relaxed"

refname = ("Danaus_R", "BMORI")
#replace refnames with name of reference for corresponding probe kit.

keep_del = []

alignIN = AlignIO.read(infile, fileformat)
outfile = "ReGAl_" + infile.split(".")[0] + ".fas"
    
for record in alignIN:
    for i in refname:
        if i in record.id:
            ref = i
            for site in record.seq:
                if site == "-":
                    keep_del.append(0)
                else:
                    keep_del.append(1)
                    
edited = alignIN[:, 0:0]
for i in range(len(keep_del)):
    if keep_del[i] == 1:
        edited += alignIN[:, i:i+1]

out_records = MultipleSeqAlignment(edited)
AlignIO.write(out_records, outfile, "fasta")