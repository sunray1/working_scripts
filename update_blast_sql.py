#!/usr/bin/env python
#script designed to take hand_changes and change them in the blast sqlite database
import os, sys, sqlite3
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
namechange = {}

with open(sys.argv[1]) as o:
    line = o.readline()
    while line:
        namechange[line.split('\t')[0].strip()] = line.split('\t')[1].strip()
        line = o.readline()
for n in namechange.keys():
    c.execute("UPDATE blast SET genus='" + namechange[n].split()[0] + "' WHERE Species='" + str(n) + "';")
    c.execute("UPDATE blast SET epithet='" + namechange[n].split()[1] + "' WHERE Species='" + str(n) + "';")
    c.execute("UPDATE blast SET Species='" + namechange[n] + "' WHERE Species ='" + str(n) + "';")

conn.commit()
conn.close()
