#! /usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ambiguous_dna, extended_protein

import numpy
import sys
import os


def align(fastain):
#	thread = arguments[arguments.index("-T")+1]
	thread = 4
	outparts=fastain.split(".")
	seqL=[]
	runfull = "mafft --adjustdirection --maxiterate 1000 --localpair --thread %s %s > %s"
	addfrag = "mafft --adjustdirection --maxiterate 1000 --localpair --thread %s --addfragments %s %s > %s"
	addfull = "mafft --adjustdirection --maxiterate 1000 --localpair --thread %s --add %s %s > %s"
	killfile = "rm %s"
	
	fullcount=0
	fragcount=0
	
	
	fasta_sequences = SeqIO.parse(open(fastain),'fasta')

	end = False
	outfile1 =outparts[0] + "_FULL.fas"
	outfile2 =outparts[0] + "_FRAGS.fas"
	outfile3 =outparts[0] + "_FULLALIGN.fa"
	outfile4 =outparts[0] + "_FINALALIGN.fa"
	
	for record in fasta_sequences:
		record.description=""
		SEQ = record.upper()
		seqL.append(len(SEQ.seq))
	
	#print seqL
	dev=numpy.std(seqL)*2
	upper=max(seqL)
	above=upper-dev
	
	f= open(outfile1, "w")
	g= open(outfile2, "w")
	
	fasta_sequences = SeqIO.parse(open(fastain),'fasta')
	end = False
	for record in fasta_sequences:
		record.description=""
		SEQ = record.upper()	
		if len(SEQ.seq) >= above:
			SeqIO.write([record], f, "fasta")		
			fullcount+=1
		else:
			SeqIO.write([record], g, "fasta")
			fragcount+=1
	
	f.close()
	g.close()
	
	print(outfile1 + ":" + str(fullcount) + " seqs\n" + outfile2 + ": " + str(fragcount) + " seqs\n")
	if fullcount == 1 and fragcount > 50:
		#if only one full length and a large amount of fragments, we want to split the fragments, align the chunks to the full length sequence.
		#might be better to align all the chunks first and then --addlong the last long sequence at the end?
		print("RUNNING 1 FULL, LOTS FRAGS")
		filelist = splitseqs(outfile2)
		for i in range(len(filelist)-1):
			prev_outfile = outparts[0] + "_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(addfrag % (thread, filelist[0], outfile1, new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))
	elif fullcount == 1 and fragcount <= 50:
		#if only one full length and not a lot of fragments, just align the fragments to the full length sequence
		print("RUNNING 1 FULL, FEW FRAGS")
		os.system(addfrag % (thread, outfile2, outfile1, outfile4))
		
	elif fullcount == 0 and fragcount > 50:
		#if there are ONLY a lot of fragments, chunk and align the fragments
		#I dont think this can happen due to how full and frags are define but in case.
		print("RUNNING 0 FULL, LOTS FRAGS")
		filelist = splitseqs(outfile2)
		for i in range(len(filelist)-1):
			print(i)
			prev_outfile = outparts[0] + "_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(runfull % (thread, filelist[0], new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))

	elif fragcount == 0 and fullcount > 50:
	#if there are ONLY a lot of full length, chunk and align the fragments
		print("RUNNING LOTS FULL, 0 FRAGS")
		filelist = splitseqs(outfile1)
		for i in range(len(filelist)-1):
			print(i)
			prev_outfile = outparts[0] + "_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(runfull % (thread, filelist[0], new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))
	elif fullcount == 0 and fragcount <= 50:
		#if there are ONLY a few fragments, align the fragments
		#I dont think this can happen due to how full and frags are define but in case.
		print("RUNNING 0 FULL, FEW FRAGS")
		os.system(runfull % (thread, outfile2, outfile4))

	elif fragcount == 0 and fullcount <= 50:
	#if there are ONLY a few full length, align the fragments
		print("RUNNING FEW FULL, 0 FRAGS")
		os.system(runfull % (thread, outfile1, outfile4))
	elif fragcount > 50 and fullcount > 50:
	#if lots of fragments and lots of full length, chunk and align both, then add the aligned frags into aligned fulls
		print("RUNNING LOTS FULL, LOTS FRAGS")
		#align fulls
		filelist = splitseqs(outfile1)
		for i in range(len(filelist)-1):
			print(i)
			prev_outfile = outparts[0] + "_FULL_aligned" + str(i) + ".fas"
			new_outfile = outparts[0] + "_FULL_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(runfull % (thread, filelist[0], new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))
		finalfullalign = new_outfile
		#align frags
		filelist = splitseqs(outfile2)
		for i in range(len(filelist)-1):
			print(i)
			prev_outfile = outparts[0] + "_FRAG_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_FRAG_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(runfull % (thread, filelist[0], new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))
		finalfragalign = new_outfile
		#align frags to full
		os.system(addfrag % (thread, finalfragalign, finalfullalign, outfile4))
	elif fragcount > 50 and fullcount <= 50:
	#if there are lots of fragments and only a few fulls (but more than 1 b/c above ifs catch first), align the full lengths, then chunk and align the fragments to the full length sequence
		print("RUNNING FEW FULL, LOTS FRAGS")
		filelist = splitseqs(outfile2)
		for i in range(len(filelist)-1):
			prev_outfile = outparts[0] + "_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_aligned" + str(i+1) + ".fas"
			if i == 0:
				if fullcount == 1:
				#if there is only one full, no need to align it.
					os.system(addfrag % (thread, filelist[i], outfile1, new_outfile))
				else:
					os.system(runfull % (thread, outfile1, outfile3))
					os.system(addfrag % (thread, filelist[i], outfile3, new_outfile))
			else:
				os.system(addfrag % (thread, filelist[i], prev_outfile, new_outfile))
	elif fullcount > 50 and fragcount <= 50:
		#if there are lots of full length seqs and only a few fragments (but more than 0)(this includes 1), chunk and align the fulls and then add the fragments to the aligned fulls
		print("RUNNING LOTS FULL, FEW FRAGS")
		filelist = splitseqs(outfile1)
		for i in range(len(filelist)-1):
			prev_outfile = outparts[0] + "_aligned" + str(i) + ".fas"
			new_outfile = outparts[0]  + "_aligned" + str(i+1) + ".fas"
			if i == 0:
				os.system(runfull % (thread, filelist[0], new_outfile))
			else:
				os.system(addfull % (thread, filelist[i], prev_outfile, new_outfile))
		i = i+1		
		prev_outfile = filelist[i].split(".")[0] + "_aligned" + i + ".fas"
		new_outfile = filelist[i].split(".")[0] + "_aligned" + i+1 + ".fas"
		os.system(addfrag % (thread, outfile2, prev_outfile, new_outfile))
	elif fullcount <= 50 and fragcount <= 50:
		print("RUNNING FEW FULL, FEW FRAGS")
		#if there are only a few of each, align the fulls and then align the fragments to the aligned full
		#and by default neither is equal to 1
		print("RUNNING MAFFT FULL ALIGNMENT.....")
		os.system(runfull % (thread, outfile1, outfile3))
		print("RUNNING MAFFT ADDING FRAGS.....")
		os.system(addfrag % (thread, outfile2, outfile3, outfile4))
	#os.system(killfile %(outfile1))
	#os.system(killfile %(outfile2))
	#os.system(killfile %(outfile3))
	
def splitseqs(infile):
    #makes and returns a list of files with only 50 sequences in them
	openfile = list(SeqIO.parse(infile, "fasta"))
	outlist = []
	flen = len(openfile)
	r = list(range(0, flen, 50))
	r.append(flen)
	count = 1
	for i in range(len(r)):
		outname = infile.split('.')[0] + "_" + str(count) + ".fa"
		outlist.append(outname)
		try:
			SeqIO.write(openfile[r[i]:r[i+1]], outname, "fasta")
			count += 1
		except:
			break
	return(outlist)
		

	
align(sys.argv[1])
