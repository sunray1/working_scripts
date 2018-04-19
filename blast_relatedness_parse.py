#!/usr/bin/python
#script takes self blast output and uses the GI numbers to look up
#in the sql blast database and make sure they're in the same family
# ./script.py blastsql taxodatabase blastout
import os, sys, sqlite3

conn = sqlite3.connect(sys.argv[1])
c = conn.cursor()
c.execute("ATTACH '" + sys.argv[2] + "' as 'db'")
hit_dic = {}
count = 0
#make a dictionary with target and query GIs from self blast
with open(sys.argv[3]) as o:
    line = o.readline()
    while line:
        hit_dic[line.split("|")[1]] = line.split("|")[5]
        line = o.readline()
#get the family for each match
#gets genus from blastsql, then uses while loop to loop up to family in taxodatabase
with open(sys.argv[3].split(".")[0] + '_errors.txt', 'w') as o:
    for i in hit_dic:
        count += 1
        count2 = 0
        count3 = 0
        query_rank = ''
        target_rank = ''
        for iter in c.execute("SELECT Species FROM blast WHERE GI='" + i + "'"):
            query_taxa = iter[0]
        #have to keep counts in case something is wrong with sql and it never find the family
        while query_rank != 'Family' and count2 < 10:
            for iter in c.execute("SELECT n2.namestr, r2.namestr FROM names n1, names_to_taxonconcepts ntt1, taxon_concepts tc1, taxon_concepts tc2, names_to_taxonconcepts ntt2, names n2, ranks r2 WHERE n1.name_id=ntt1.name_id AND ntt1.tc_id=tc1.tc_id AND tc1.parent_id=tc2.tc_id AND ntt2.tc_id=tc2.tc_id AND n2.name_id=ntt2.name_id AND tc2.rank_id = r2.rank_id AND n1.namestr='"+query_taxa+"' GROUP BY n2.namestr, r2.namestr"):
                query_taxa = str(iter[0])
                query_rank = str(iter[1])
                count2 += 1
        if count2 == 10:
            o.write('ERROR WITH ' + str(i) + '/' + str(hit_dic[i])+'\n')
        for iter in c.execute("SELECT Species FROM blast WHERE GI='" + hit_dic[i] + "'"):
            target_taxa = iter[0]
        while target_rank != 'Family' and count3 < 10:
            for iter in c.execute("SELECT n2.namestr, r2.namestr FROM names n1, names_to_taxonconcepts ntt1, taxon_concepts tc1, taxon_concepts tc2, names_to_taxonconcepts ntt2, names n2, ranks r2 WHERE n1.name_id=ntt1.name_id AND ntt1.tc_id=tc1.tc_id AND tc1.parent_id=tc2.tc_id AND ntt2.tc_id=tc2.tc_id AND n2.name_id=ntt2.name_id AND tc2.rank_id = r2.rank_id AND n1.namestr='"+target_taxa+"' GROUP BY n2.namestr, r2.namestr"):
                target_taxa = str(iter[0])
                target_rank = str(iter[1])
                count3 += 1
        if count3 == 10:
            o.write('ERROR WITH ' + str(hit_dic[i]) + '/' + str(i)+'\n')
        print(round((float(count)/float(len(hit_dic.keys())))*100, 2))
    #if the families are not the same between the query and the hit, print
        if query_taxa != target_taxa and count2 != 10 and count3 != 10:
            o.write(i+' '+query_taxa+' '+hit_dic[i]+' '+target_taxa+'\n')
conn.close()