#!/usr/bin/env python


#./rename_seqs fasta blast_results.db tax.db Order_Family_Species

from Bio import SeqIO
import sys, sqlite3


filein = sys.argv[1]
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
c.execute("ATTACH '" + sys.argv[3] + "' as 'tax'")
template = sys.argv[4].split("_")

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
    try:
        accnum = accnum1.split("_R_")[-1]
    except:
        accnum = accnum1
    for iter in c.execute("SELECT tc_id, Species FROM blast WHERE accession = '" + str(accnum) + "'"):
        tc_id= iter[0]
        species = str(iter[1]).replace(" ", "_")
    rank = ''
    taxonomy = []
    ranks = []
    count1 = 0
    sp_tc_id = tc_id
    while rank != 'Order' and count1 != 20:
        for iter in c.execute("SELECT r.namestr, tc.parent_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "' AND tc.rank_id = r.rank_id"):
            sp_tc_id = iter[1]
            rank = iter[0]
            ranks.append(rank)
        for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "' and validity = 'valid'"):
            taxonomy.append(str(iter[0]))
        count1 += 1
    if count1 == 20:
        pass
    else:          
        for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(tc_id) + "' and validity = 'valid'"):
            species = iter[0]
        taxonomy = [species.replace(' ', '_')] + taxonomy[:len(taxonomy)-1]
        new_name = []
        for i in template:
            if i in ranks:
                new_name.append(taxonomy[ranks.index(i)])
            else:
                new_name.append('X')
        
        
        record.id = '_'.join(new_name)
    #print(record.id)
    record.description = ''
    recs_changed.append(record)
outfile = open(filein.split(".")[0] + "_renamed.fa", "w")
SeqIO.write(recs_changed, outfile, "fasta")
