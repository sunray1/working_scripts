#!/usr/bin/env python
#use this when sequences are already aligned
from Bio import SeqIO, AlignIO
from Bio import SeqIO, AlignIO
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import sys
muscle_cline = MuscleCommandline(clwstrict=True)
filein = sys.argv[1]

seq_dic = {}
seqOUT = []
outfile = filein.rsplit(".", 1)[0] + "_dupsmerge.fas"

records = AlignIO.read(filein, "fasta")
for record in records:
    try:
        seq_dic[record.id].append(record)
    except:
        seq_dic[record.id] = [record]
        
for i in seq_dic:
    if len(seq_dic[i]) > 1:
        align = MultipleSeqAlignment(seq_dic[i])
        # handle_string = StringIO()
        # SeqIO.write(seqs, handle_string, "fasta")
        # data = handle_string.getvalue()
        # stdout, stderr = muscle_cline(stdin=data)
        # align = AlignIO.read(StringIO(stdout), "clustal")

                
        count = 0
        gaps = 0
        for col in range(0, len(align[0])):
            column = align[:,col].replace("N", "").replace("-", "")
            if len(column) > 1:
                if column[1:]==column[:-1]:
                    count=count+1
            else:
                gaps=gaps+1
            
        if len(align[0])-gaps == 0:
            iden = 0.0
        else:
            iden = 100*(count/float((len(align[0])-gaps)))
        print(iden)
        
        #change Ns with data in other columns to gaps so consensus will work
        data = ["A", "C", "T", "G"]
        for x, record in enumerate(align):
            record.seq = record.seq.tomutable()
            for col in range(0, len(align[0])):
                column = align[:,col]
                if column[x] == "N":
                    if any(nt in column[:x] + column[x+1:] for nt in data) == True:
                        record.seq[col] = "-"
                        
        if iden > 90 or iden == 0:
            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.dumb_consensus(ambiguous="N", threshold = .49)
            seq_rec = SeqRecord(consensus)
            seq_rec.id = i
            seq_rec.description = ""
            seqOUT.append(seq_rec)
        else:
            seqs = seq_dic[i]
            for rec in seqs:
                rec.description = ""
                seqOUT.append(rec)
    else:
        seqs = seq_dic[i]
        for rec in seqs:
            rec.description = ""
            seqOUT.append(rec)

SeqIO.write(seqOUT, outfile, "fasta")
