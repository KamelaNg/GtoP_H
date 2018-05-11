#!/usr/bin/env python
from os import system 
import glob 
import sys
import re
from collections import defaultdict
from sets import Set
from itertools import groupby
import argparse
import copy

"""
Author: Kamela Ng <11 May 2018>

This script converts whole genome sequence data in the form of tab files into the most liekly output that would be observed from a GenoType MTBDRplus v2.0  (hereinafter called Hainplus) run on this sample

Input:
A folder containing all the vcf files that are to be processed
A map file that states the relationship between codons where mutations are known to be detected by Nipro (based on this following paper: Ng KC, Meehan CJ, Torrea G, Goeminne L, Diels M, Rigouts L, de Jong BC, Andre E
Potential Application of Digitally Linked Tuberculosis Diagnostics for Real-Time Surveillance of Drug-Resistant Tuberculosis Transmission: Validation and Analysis of Test Results. JMIR Med Inform 2018;6(1):e12.) and positions in the genome
A filename for the output

The map file is optional. If not supplied, it is assumed the standard H37Rv NC000962.3 was used and that the file sample_mapfile.txt is in the current working folder
An example of this file is (tab separated columns):
Codon	Ref	Gen	Pos	
428	761088
430	761094
431	761097
432	761100
434	761106
435	761109
437	761115
441	761127
445	761139
446	761142
450	761154
452	761160

The filename for the output is optional. If not supplied, output will be genomeToProbe_Hainplus_2.txt

Output:
The output file contains, for each input vcf file:
resistance or susceptibility to rifampicin labelled as RIF resistance detected or not detected, mutant codon position, mutant codon nucleotides,  absent wildtype and developing mutant probes, and the reactions for each probe corresponding the specific mutation in the sample - with 1 representing developing probe and 0 representing absent probe.

Usage:
python Hainplus_GtoP_final.py --folder <folderName> --map <codonMapFile> (optional) --out <outputFilename> (optional)
"""
#Hardcode the codon to wildtype nucleotides
codon_nuc = {
			'428' : ['A', 'G', 'C'], 
			'430' : ['C', 'T', 'G'], 
			'431' : ['A', 'G', 'C'], 
			'432' : ['C', 'A', 'A'], 
			'434' : ['A', 'T', 'G'], 
			'435' : ['G', 'A', 'C'], 
			'437' : ['A', 'A', 'C'], 
			'441' : ['T', 'C', 'G'], 
			'445' : ['C', 'A', 'C'], 
			'446' : ['A', 'A', 'G'],
			'450' : ['T', 'C', 'G'], 
			'452' : ['C', 'T', 'G']
}

#hardcode the codon positions with the mutated codons and the associated capturing probe
Hainplusmut = {
			'428' : {'AGG':[['WT1'],['']]},
			'430' : {'CCG':[['WT2'],['']]}, 
			'431' : {'GGC':[['WT2'],['']]}, 
			'432' : {'GAA':[['WT2','WT3','WT4'],['']]},
			'434' : {'GTG':[['WT3'],['']],'ATT':[['WT3'],['']],'ACG':[['WT3'],['']]}, 
			'435' : {'GTC':[['WT3','WT4'],['MUT1']],'TAC':[['WT3','WT4'],['']],'TTC':[['WT3','WT4'],['']],'GAA':[['WT3','WT4'],['']]}, 
			'437' : {'GAC':[['WT4'],['']]},
			'441' : {'CAG':[['WT5','WT6'],['']],'TTG':[['WT5','WT6'],['']]},
			'445' : {'TAC':[['WT7'],['MUT2A']],'GAC':[['WT7'],['MUT2B']],'GGC':[['WT7'],['']],'AGC':[['WT7'],['']],'TCC':[['WT7'],['']],'AAC':[['WT7'],['']],'CAG':[['WT7'],['']],'CAA':[['WT7'],['']]}, 
			'446' : {'CAG':[['WT7'],['']]},
			'450' : {'TTG':[['WT8'],['MUT3']],'TGG':[['WT8'],['']],'TTC':[['WT8'],['']]}, 
			'452' : {'CCG':[['WT8'],['']]}
}

			
#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True, help='Folder containing the tab files to be processed')
parser.add_argument('--map', required=False, default="sample_mapfile.txt", help='File listing the relationship between codon and genome position (default is "sample_maplfile.txt")')
parser.add_argument('--out', required=False, default="genomeToProbe_Hainplus_2.txt", help='Filename for output (default is "genomeToProbe_Hainplus_2.txt")')

args = parser.parse_args()
	
#read in map file 
try:
	mapf=open(args.map, 'rU')
except IOError:
	print "\n map file not found."
	sys.exit()

#Open output files
try:
	HainplusOutputfinal = open(args.out, "w")  
except IOError:
	print 'no room for save file'						
	sys.exit()
HainplusOutputfinal.write("Filename"+"\t"+"RIF Resistance"+"\t"+"MutCodonNum"+"\t"+"MutCodonVal"+"\t"+"AbsentProbe"+"\t"+"DevelopingProbe"+"\t"+"WT1"+"\t"+"WT2"+"\t"+"WT3"+"\t"+"WT4"+"\t"+"WT5"+"\t"+"WT6"+"\t"+"WT7"+"\t"+"WT8"+"\t"+"MUT1"+"\t"+"MUT2A"+"\t"+"MUT2B"+"\t"+"MUT3"+"\n")
  
#read in map and create a dictionary of the codons to genome positions
#each key is the three codon positions (supplied genome +0, +1, +2) and the value is the codon
mapdict={}
while 1:
	line=mapf.readline()
	if not line:
		break
	line=line.rstrip()
	if re.match("Codon",line):#on an info line so skip
		continue
	sections=line.split("\t")
	pos=sections[1]
	mapdict[pos, str(int(pos)+1), str(int(pos)+2)]=sections[0]
	codonpos=mapdict.keys() #get a list of all the genome positions, grouped as triplets
mapf.close()

#output probe pattern
outputList=["1","1","1","1","1","1","1","1","0","0","0","0"]
outputKeys=["WT1","WT2","WT3","WT4","WT5","WT6","WT7","WT8","MUT1","MUT2A","MUT2B","MUT3"]

#read in vcf file
#for each file make a copy of the WT codon patterns, mutate it based on the vcf file and then output the associated result
OpenDir = glob.glob(args.folder + "/*")										  
for File in OpenDir: 
	#print File 
	if File.endswith(".vcf"): 
		#print File
		try:
			vcf=open(File,'rU')
			#print tab
		except IOError:
			print "\n VCF file not found."
		sampleName=File.rsplit(".",1)[0].split("/")[-1]
		print 'Processing '+sampleName
		
		sampleCodons=copy.deepcopy(codon_nuc)
		while 1:
			line=vcf.readline()
			#print s
			if not line:
				break
			line=line.rstrip()
			#print s
			if re.match("#",line):	#on an info line so skip
				continue
			sections=line.split("\t")
			#print sections
			mutpos=sections[1]
			altbase=sections[4]
			#print mutpos, altbase
			#go through the genome positions that are associated with any change in a codon base
			for sublist in codonpos:
				if mutpos in sublist: #check if the genome position from the vcf file is a position in a codon of interest
					codon=mapdict[sublist]
					#print codon
					if codon in Hainplusmut.keys(): #check if the codon is in the mutation list
						Hainplusmutcodval = sampleCodons[codon]
						ind = sublist.index(mutpos)
						#modify the WT codon to be the new mutated codon
						Hainplusmutcodval[ind] = altbase
						sampleCodons[codon]=Hainplusmutcodval					
		vcf.close()
			
		#output results
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		wtprobeList=[]
		mutprobeList=[]
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		codonList=codon_nuc.keys()
		codonList.sort()
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos]: #the codons are not the same 
				codonCheck="".join(sampleCodons[codonPos])
				codnucres=Hainplusmut[codonPos]
				#print codnucres
				if codonCheck in codnucres:
					mutCodonValList.append(codonCheck)
					mutCodonList.append(codonPos)
					for wtprobe in (codnucres[codonCheck][0]):
						sampleOutputList[outputKeys.index(wtprobe)]="0"
						if wtprobe not in wtprobeList:
							wtprobeList.append(wtprobe)
					for mutprobe in (codnucres[codonCheck][1]):	
						if mutprobe in outputKeys:
							sampleOutputList[outputKeys.index(mutprobe)]="1"						
							mutprobeList.append(mutprobe)
						
							
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			HainplusOutputfinal.write(sampleName+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(wtprobeList)+"\t"+",".join(mutprobeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			HainplusOutputfinal.write(sampleName+"\t"+"NOT DETECTED\n")
HainplusOutputfinal.close()
sys.exit()