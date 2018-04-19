#!/usr/bin/env python

from Bio import Entrez, SeqIO
import os, sys

Entrez.email = "sunray1@ufl.edu"
seqids = []
records = []

if len(sys.argv) < 2:
	sys.exit("module load python/2.7.6\ngetseqsfromGenbank.py id_list")
name = sys.argv[1].split(".")[0]
with open(sys.argv[1]) as o:
	line = o.readline()
	while line:
		seqids.append(line.strip())
		line = o.readline()
seqlists = [seqids[i:i+200] for i in range(0, len(seqids), 200)]
for i in seqlists:
	seqids_sub = ",".join(i)
	handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=seqids_sub)
	for seq_record in SeqIO.parse(handle, "fasta"):
		records.append(seq_record)
	print(str(round((float(len(records))/float(len(seqids)))*100, 2)) + "%")
	handle.close()

SeqIO.write(records, name + ".fa", "fasta")

