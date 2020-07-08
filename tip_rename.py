#!/usr/bin/env python

import sys, os
assert sys.version_info >= (3, 0), "Needs python3 or higher"
from io import StringIO #python3
try:
    from Bio import Phylo
except:
    sys.exit("Error: Biopython needs to be installed. - try module load python3")
from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd



argp = ArgumentParser(description='Renames tip labels on a phylogenetic tree', formatter_class=RawTextHelpFormatter)
argp.add_argument('-i', '--input', help='Input tree', required=True)
argp.add_argument('-t', '--taxonomy', help='Adds higher taxonomy to tips - csv taxonomy file')
argp.add_argument('-tf', '--taxformat', help='List of what columns from taxonomy to add (will add before species name), ie. "Family, Subfamily, Tribe" (default)')
argp.add_argument('-l', '--loci', help='Adds number of loci to tips - csv loci number file')
argp.add_argument('-s', '--subset', help='Only change a subset of the tips - if no file, will rename all tips')
argp.add_argument('-r', '--replace', help="Change tips based on a before/after file - csv file")
argp.add_argument('-f', '--format', help='Format of tree file', choices={"newick", "nexus"}, required=True)
argp.add_argument('-o', '--output', help='Output tree, will add \'_renamed\' as prefix if not specified')
argp.add_argument('-n', '--names', help='Prints an outfile with the labels of the input tree', action="store_true")
args = argp.parse_args()
tformat = args.format
input_tree = args.input

if len(sys.argv) == 1:
    argp.print_help()
    sys.exit()
    

    
tree = Phylo.read(input_tree, tformat)

#to print tip names
if args.names:
    if not args.output:
        output = input_tree.split(".")[0] + "_tipnames.txt"
    else:
        output = args.output
    term_names = [term.name for term in tree.get_terminals()]
    with open(output, "w") as o:
        for i in term_names:
            o.write(i + "\n")
outstrand = StringIO()
Phylo.write(tree, outstrand, "newick")
newick_tree = outstrand.getvalue()

if args.subset:
    print("Only changing a subset of the taxonomy found in " + args.subset)
    tip_names = [line.strip() for line in open(args.subset)]
else:
    tip_names = [tip.name for tip in tree.get_terminals()]
    
#parse taxonomy    
if args.taxonomy:
    print("Parsing taxonomy...")

    if not args.taxformat:
        argp.print_help()
        sys.exit("Need taxonomy format - ie \"Family, Subfamily, Tribe\"")    
    tax_dict = {}
    synonym_dict = {}
    tax_file = args.taxonomy
    tax_df = pd.read_csv(tax_file, header=0)
    taxformat = args.taxformat.split(", ")
    #iterate through each row
    for row in tax_df.iterrows():
        row_headers = []
        #parse out taxonomy
        for header in taxformat:
            if type(row[1][header]) == str:
                row_headers.append(row[1][header].title())
            row_taxonomy = "_".join(row_headers)
        #add species to dictionary
        species = row[1]["Species"].replace(" ", "_")
        tax_dict[species] = row_taxonomy + "_" + species
        #add synonyms to dictionary
        if type(row[1]["Synonyms"]) == str:
            for syn in row[1]["Synonyms"].split(", "):
                syn_key = syn.split()[0]
                synonym_dict[syn_key] = row_taxonomy + "_" + syn_key
    print("Adding higher taxonomy...")
    tip_str = ",".join(tip_names)
    
    for key in tax_dict.keys():
        tip_str = tip_str.replace(key, tax_dict[key])

    tip_names_replaced = tip_str.split(",")
    for num, tip_replace in enumerate(tip_names_replaced):
        newick_tree = newick_tree.replace(tip_names[num], tip_replace)

if args.loci:
    print("Adding numbers of loci found in " + args.loci)
    loci_dict = {}
    loci_name_list = []
    loci_num_list = []
    loci_file = args.loci
    loci_df = pd.read_csv(loci_file)
    #iterate through each row
    for row in loci_df.iterrows():
        loci_name_list.append(row[1][0])
        loci_num_list.append(row[1][1])
    
    #use higher taxonomy if there is a taxonomy
    
    if args.taxonomy:
        loci_name_str = ",".join(loci_name_list)
        for key in tax_dict.keys():
            loci_name_str = loci_name_str.replace(key, tax_dict[key])
        loci_name_replaced = loci_name_str.split(",")
        for n, i in enumerate(loci_name_replaced):
            loci_dict[i] = i + "_" + str(loci_num_list[n])
    
    if not args.taxonomy:
        for n, i in enumerate(loci_name_list):
            loci_dict[i] = i + "|" + str(loci_num_list[n])
    
    for loci in loci_dict:
        newick_tree = newick_tree.replace(loci, loci_dict[loci])
    
    if args.output:
        output = args.output
    else:
        output = input_tree.split(".")[0] + "_renamed." + input_tree.split(".")[1]
    
    with open(output, "w") as out:
        out.write(newick_tree)
        
if args.replace:
    print("Replacing tip names found in " + args.replace)
    replace_df = pd.read_csv(args.replace, header=0)
    #iterate through each row
    for row in replace_df.iterrows():
        tip_from = row[1][0]
        tip_to = row[1][1]
        newick_tree = newick_tree.replace(tip_from, tip_to)
        
    if args.output:
        output = args.output
    else:
        output = input_tree.split(".")[0] + "_renamed." + input_tree.split(".")[1]
    
    with open(output, "w") as out:
        out.write(newick_tree)
        
        