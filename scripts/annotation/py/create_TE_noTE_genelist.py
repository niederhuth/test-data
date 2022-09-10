#! /usr/bin/python
#Original by Tiffany Liu 11/17/2016
#Modified by Chad Niederhuth 6/1/2022

import argparse
import os.path
import sys

#Add and parse command line arguments
parser = argparse.ArgumentParser(description="This filters an input gene list for putative TE genes")
parser.add_argument("--input_geneList", help="Input list of genes, where each line has a gene name.", required = True)
parser.add_argument("--pfamhmm", help="hmmscan output of pfam domains (requires --TEpfam_list", required = False)
parser.add_argument("--TEpfam_list", help="List of TE-related Pfam IDs.", required = False)
parser.add_argument("--TEhmm", help="hmmscan output of TE-related genes (e.g. gypsy hmm)", required = False)
parser.add_argument("--TEblast", help="Blast results against TE sequences (requires outfmt 6).", required = False)
parser.add_argument("--output_non_TE_genes", help="Path to output file for non TE genes.", required = True)
parser.add_argument("--output_TE_like_genes", help="Path to output file for possible TE-like genes.", required = True)
args = parser.parse_args()

#Check for input_geneList
input_geneList = args.input_geneList
if os.path.isfile(input_geneList) == False:
	print('ERROR: he file ' + input_geneList + ' does not exist.\n')
	sys.exit()

#Create an empty set of TE genes to be filtered
TEgeneSet = set()

#Filter based on pfam hits
if args.pfamhmm:
	#Check to make sure Tpfam_list is set
	if not args.TEpfam_list:
		print('ERROR: --pfamhmm requires --TEhpfam_list')
		sys.exit()
	#Check for pfamhmm input
	pfamhmm = args.pfamhmm
	if os.path.isfile(pfamhmm) == False:
		print('ERROR: The file ' + pfamhmm + ' does not exist.\n')
		sys.exit()
	#Check for TEpfam_list file
	TEpfam_list = args.TEpfam_list
	if os.path.isfile(TEpfam_list) == False:
		print('ERROR: The file ' + TEpfam_list + ' does not exist.\n')
		sys.exit()
	#Read in list of TE-related pfam domains
	TEpfamSet = set()
	with open(TEpfam_list) as input_fh_TEpfam:
		for each_line in input_fh_TEpfam:
			if 'Pfam' == each_line[0:4]:
				pass #Ignoring header
			else:
				line_string = each_line.strip()
				(Pfam,Domain,Description) = line_string.split('\t')
				TEpfamSet.add(Pfam)
	#Compare PfamIDs between the TEpfam_list and the pfamhmm results and remove genes with TE pfam domains
	pfamSet = set()
	with open(pfamhmm) as input_fh_maxPfam:
		for each_line in input_fh_maxPfam:
			if "#" not in each_line[0]: #ignoring header
				line_string = each_line.strip()
				line_tuple = line_string.split()
				pfam = line_tuple[1]
				pfamSplit = pfam.split(".")
				pfamID = pfamSplit[0]
				geneID = line_tuple[2]
				if pfamID in TEpfamSet:
					TEgeneSet.add(geneID)
					pfamSet.add(geneID)
	output_pfam = open("pfam_filtered_genes.txt", 'w') 
	for each_element in pfamSet:
		output_pfam.write("%s\n"%(each_element))
	output_pfam.close()    
	print("Number of TE-related genes identified by pfam hmmscan: ", len(pfamSet))


#Filter based on TE hmm results
if args.TEhmm:
	#Check for TEhmm file
	TEhmm = args.TEhmm
	if os.path.isfile(TEhmm) == False:
		print('ERROR: The file ' + TEhmm + ' does not exist.\n')
		sys.exit()
	#Filter genes identified by TE hmmscan
	TEhmmSet = set()
	with open(TEhmm) as input_fh_TEhmm:
		for each_line in input_fh_TEhmm:
			if "#" not in each_line[0]:
				line_string = each_line.strip()
				line_tuple = line_string.split()
				gene = line_tuple[2]
				TEgeneSet.add(gene)
				TEhmmSet.add(gene)
	output_TEhmm = open("TEhmm_filtered_genes.txt", 'w') 
	for each_element in TEhmmSet:
		output_TEhmm.write("%s\n"%(each_element))
	output_TEhmm.close()    
	print("Number of TE-related genes identified by TE hmmscan: ", len(TEhmmSet))            

#Filter based on TE blast results
if args.TEblast:
	#Check for TEblast file
	TEblast = args.TEblast
	if os.path.isfile(TEblast) == False:
		print('\n\nThe file ' + TEblast + ' does not exist.\n')
		sys.exit()
	#Filter genes identified by TE hmmscan
	blastSet = set()
	with open(TEblast) as input_fh_TEblast:
		for each_line in input_fh_TEblast:
			if "#" not in each_line[0]:
				line_string = each_line.strip()
				line_tuple = line_string.split()
				gene = line_tuple[0]
				TEgeneSet.add(gene)
				blastSet.add(gene)
	output_blast = open("blast_filtered_genes.txt", 'w') 
	for each_element in blastSet:
		output_blast.write("%s\n"%(each_element))
	output_blast.close()    
	print("Number of TE-related genes identified by blast: ", len(blastSet)) 

print("Total number of TE-related genes identified: ", len(TEgeneSet))

#Add input-geneList as a set
oldGeneSet = set()
with open(input_geneList) as input_fh_geneList:
	for each_line in input_fh_geneList:
		line_string = each_line.strip()
		oldGeneSet.add(line_string)

#Remove genes TE-related genes from gene list and output new gene list
newGeneSet = oldGeneSet.difference(TEgeneSet)
print("Length of original gene list: ", len(oldGeneSet))
print("Length of non-TE gene list: ", len(newGeneSet))
#Output list of non-TE genes
output_non_TEs_fh = open(args.output_non_TE_genes, 'w')
for each_element in newGeneSet:
	output_non_TEs_fh.write("%s\n"%(each_element))
output_non_TEs_fh.close()
#Output list of putative TE genes
output_TEs_fh = open(args.output_TE_like_genes, 'w')
for each_element in TEgeneSet:
	output_TEs_fh.write("%s\n"%(each_element))
output_TEs_fh.close()








