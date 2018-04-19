#!/usr/bin/env python
#script designed to check wikipedia for synonyms, searches butterflynet database for matching validnames
import os, sys, urllib2, sqlite3
from ast import literal_eval
names = []
user_agent = 'Mozilla/5.0 (Windows NT 6.1; Win64; x64)'
API_names = "https://api.mol.org/0.x/taxonomy/tcsearch?searchstrs=%s"
API_tcid = "https://api.mol.org/0.x/taxonomy/tcinfo?tcids=%s&getsynonyms=True"
wiki = "https://en.wikipedia.org/wiki/%s"
headers = {'User-Agent': user_agent}
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
#c.execute("ATTACH '" + sys.argv[2] + "' as 'db'")
for line in open(sys.argv[1]):
	names.append(line.strip())
with open("API_edit.txt", "w") as o:
	for i in names:
		names2 = set()
		ids = set()
		syns = []
		name = i.replace(" ", "_")
		try:
			print(wiki % (name))
			req = urllib2.Request(wiki % (name),None,headers)
			response = urllib2.urlopen(req)
			data = response.read()
	#		o.write(data)
			spl1 = data.split('Synonyms')[1]
			spl2 = spl1.split('</table>')[0]
			spl3 = spl2.split('<i>')[1:]
			for l in spl3:
				syns.append(l.split('</i>')[0])
			for n in syns:
				for iter in c.execute("SELECT ntt.tc_id FROM names n, names_to_taxonconcepts ntt WHERE n.namestr = '" + n + "' and n.name_id = ntt.name_id;"):
					ids.add(iter)
				for r in ids:
					for iter2 in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE tc_id = " + str(iter).replace("(", "").replace(",)", "") + " and n.name_id = ntt.name_id and ntt.validity = 'valid';"):
						names2.add(str(iter2))
			o.write(i + '\t' + str(list(names2)) + '\n')
#		if len(names2) > 0:
#			o.write(i + '\t' + str(n) + '\n')
		except:
			pass
#		dic = literal_eval(data)
#		print(dic.keys())
#		for i in dic.values()[0]:
#			req = urllib2.Request(API_tcid % (i['tcid']),None,headers)
#			response = urllib2.urlopen(req)
#			data = response.read()
#			tcid_dic = literal_eval(data)
#			o.write(str(dic.keys()) + '\t' + str(tcid_dic.values()[0]['name']) + '\n')
#			print('syn' + str(tcid_dic.values()[0]['synonyms']))
conn.close()