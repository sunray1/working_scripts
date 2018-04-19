#!/usr/bin/env python
# This script takes the matrix that comes out of a pairwise comparison and describes the
# percentage of species in each family that has 3 or more genes
import os, sys, sqlite3
family_list = []
with open(sys.argv[1]) as o:
    line = o.readline()
    species_list = line.replace(".", " ").replace("\"", "")
species_list =species_list.strip().split(',')[1:]
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
#first get a list of all families
for iter in c.execute("SELECT ntt.tc_id FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE r.namestr = \"Family\" and n.name_id=ntt.name_id and tc.tc_id=ntt.tc_id and tc.rank_id=r.rank_id"):
    family_list.append(iter[0])
for fam in family_list:
    idnums=[fam]
    typelist = ["Family"]
    gen_dic = {}
    sp_in_genus = []
    #for each family, get all the species
    for n in idnums:
        if typelist[idnums.index(n)] != 'Species':
                for iter in c.execute("SELECT tc.tc_id, r.namestr FROM taxon_concepts tc, ranks r WHERE tc.rank_id=r.rank_id AND tc.parent_id =" + str(n) + ";"):
                    idnums.append(iter[0])
                    typelist.append(str(iter[1]))
    id_type_dic = dict(zip(idnums, typelist))
# pull out everything that aren't species and pull out genera - these are used to pull out species from the matrix
    for n in id_type_dic.keys():
        if id_type_dic[n] == "Genus":
            for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE tc.tc_id = "+ str(n) + " and n.name_id=ntt.name_id and tc.tc_id=ntt.tc_id and tc.rank_id=r.rank_id and ntt.validity = 'valid'"):
               gen_dic[str(iter[0])] = "Genus"
        if id_type_dic[n] != "Species":
           del id_type_dic[n]
#pulling out species from the matrix that is in that family
    for x in species_list:
        if x.split()[0] in gen_dic.keys():
            sp_in_genus.append(x)
    for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE tc.tc_id = "+ str(fam) + " and n.name_id=ntt.name_id and tc.tc_id=ntt.tc_id and tc.rank_id=r.rank_id and ntt.validity = 'valid'"):
        fam_str = str(iter[0])
    print("Number of Species in Family " + fam_str + " with 3 or more genes: " + str(len(sp_in_genus)))
    print("Total number of Species in Family " + fam_str + ": " + str(len(id_type_dic.keys())))
    print("% of Species within Family that have 3 or more: " + str(float(len(sp_in_genus))/float(len(id_type_dic.keys()))))
    print("% of Species within those with 3 or more in "+ fam_str+": " + str(float(len(sp_in_genus))/float(len(species_list))))
    print("Accounting for database size: " + str(((float(19242)/float(7))/float(len(id_type_dic.keys())))*float(len(sp_in_genus))/float(2594))+"\n")
conn.close()
