#!/usr/bin/env python

from Bio import SeqIO
import sys, sqlite3


filein = sys.argv[1]
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
c.execute("ATTACH '" + sys.argv[3] + "' as 'tax'")
recs_changed = []
count = 1
records = list(SeqIO.parse(filein, "fasta"))
for record in records:
    print(float(count)/float(len(records)))
    count += 1
    accnum1 = record.id.split(' ')[0]
    try:
        accnum = accnum1.split("/")[0]
    except:
        accnum = accnum1
    for iter in c.execute("SELECT tc_id, Species FROM blast WHERE accession = '" + str(accnum) + "'"):
        tc_id= iter[0]
        species = str(iter[1]).replace(" ", "_")
    rank = ''
    taxonomy = []
    sp_tc_id = tc_id
    while rank != 'Family':
        for iter in c.execute("SELECT r.namestr, tc.parent_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "' AND tc.rank_id = r.rank_id"):
            sp_tc_id = iter[1]
            rank = iter[0]
        for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "' and validity = 'valid'"):
            taxonomy.append(str(iter[0]))
    family = taxonomy[-2]
    for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(tc_id) + "' and validity = 'valid'"):
        species = iter[0]
    record.id = family + "_" + species.replace(" ", "_")
    #print(record.id)
    record.description = ''
    recs_changed.append(record)
outfile = open(filein.split(".")[0] + "_renamed.fa", "w")
SeqIO.write(recs_changed, outfile, "fasta")