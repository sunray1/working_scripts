#!/usr/bin/env python

import os
from argparse import ArgumentParser
from cleanlib.getseqsfromGenbank import pullseqs
from alignlib.process_files import process
from alignlib.rename_seqs import rename
from alignlib.add_dup_nums import number
from alignlib.LONGREF_ALIGN import align

argp = ArgumentParser(description='Takes all files, pulls down sequences and cats them together, then renames and aligns')
argp.add_argument('-f', '--files', help='list of files to be processed', nargs="*")
argp.add_argument('-b', '--blastdb', help='blastdb used to split GIs')
argp.add_argument('-g', '--fragmentgene', help='name of gene that may be fragmented')
argp.add_argument('-t', '--taxdb', help='taxonomy database')
args = argp.parse_args()
files = args.files

print('Processing files')
process(files, args.blastdb, pullseqs)

files_to_rename = [f for f in os.listdir(".") if f.endswith("_GIs.fa")]
files_to_rename.remove(str(args.fragmentgene) + "_GIs.fa")
files_to_rename.append(str(args.fragmentgene) + "_GIs_all.fa")
for i in files_to_rename:
    print("Renaming " + str(i))
    rename(i, args.blastdb, args.taxdb)
    
files_to_number = [f for f in os.listdir(".") if f.endswith("_renamed.fa")]
for i in files_to_number:
    print("Numbering " + str(i))
    number(i)
    
files_to_align = [f for f in os.listdir(".") if f.endswith("_renamed_with_nums.fa")]
for i in files_to_align:
    print("Aligning " + str(i))
    align(i)