#!/usr/bin/env python
import os, sys, sqlite3
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
count = 1
dic = {}
with open(sys.argv[1]) as o:
    line = o.readline()
    while line:
        print(float(count)/float(11868))
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI='"+ line.strip() + "'"):
            gene = str(iter[0])
        if gene in dic.keys():
            dic_list = dic[gene]
            dic_list.append(line.strip())
            dic[gene] = dic_list
        else:
            dic[gene] = [line.strip()]
        line = o.readline()
        count += 1
for g in dic:
    with open(g+"_GIs.txt", "w") as o:
        for m in dic[g]:
            o.write(m+'\n')
        