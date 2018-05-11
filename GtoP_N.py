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
Author: Kamela Ng <date>

This script converts whole genome sequence data in the form of tab files into the most liekly output that would be observed from a GenoscholarNTM+MDRTBII  (hereinafter called Nipro) run on this sample

Input:
A folder containing all the vcf files that are to be processed
A map file that states the relationship between codons where mutations are known to be detected by Nipro (based on this following paper: Ng KC, Meehan CJ, Torrea G, Goeminne L, Diels M, Rigouts L, de Jong BC, Andre E
Potential Application of Digitally Linked Tuberculosis Diagnostics for Real-Time Surveillance of Drug-Resistant Tuberculosis Transmission: Validation and Analysis of Test Results
JMIR Med Inform 2018;6(1):e12
URL: http://medinform.jmir.org/2018/1/e12/
doi:10.2196/medinform.9309
PMID: 29487047) and positions in the genome
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

The filename for the output is optional. If not supplied, output will be genomeToProbe_Nipro.txt

Output:
The output file contains, for each input tab file:
xxxxx

Usage:
python Nipro_GtoP_final.py --folder <folderName> --map <codonMapFile> (optional) --out <outputFilename> (optional)
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
Nipromut = {
			'428' : {'AGG':[['S1'],['']]},
			'430' : {'CCG':[['S1'],['']]}, 
			'431' : {'GGC':[['S1'],['']]}, 
			'432' : {'GAA':[['S1'],['']]},
			'434' : {'GTG':[['S2'],['']],'ATT':[['S2'],['']],'ACG':[['S2'],['']]}, 
			'435' : {'GTC':[['S2'],['R2']],'TAC':[['S2'],['']],'TTC':[['S2'],['']],'GAA':[['S2'],['']]}, 
			'437' : {'GAC':[['S2'],['']]},
			'441' : {'CAG':[['S3'],['']],'TTG':[['S3'],['']]},
			'445' : {'TAC':[['S4'],['R4a']],'GAC':[['S4'],['R4b']],'GGC':[['S4'],['']],'AGC':[['S4'],['']],'TCC':[['S4'],['']],'AAC':[['S4'],['']],'CAG':[['S4'],['']],'CAA':[['S4'],['']]}, 
			'446' : {'CAG':[['S4'],['']]},
			'450' : {'TTG':[['S5'],['R5']],'TGG':[['S5'],['']],'TTC':[['S5'],['']]}, 
			'452' : {'CCG':[['S5'],['']]}
}

			
#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True, help='Folder containing the tab files to be processed')
parser.add_argument('--map', required=False, default="sample_mapfile.txt", help='File listing the relationship between codon and genome position (default is "sample_maplfile.txt")')
parser.add_argument('--out', required=False, default="genomeToProbe_Nipro.txt", help='Filename for output (default is "genomeToProbe_Nipro.txt")')

args = parser.parse_args()
	
#read in map file 
try:
	mapf=open(args.map, 'rU')
except IOError:
	print "\n map file not found."
	sys.exit()

#Open output files
try:
	NiproOutputfinal = open(args.out, "w")  
except IOError:
	print 'no room for save file'						
	sys.exit()
NiproOutputfinal.write("Filename"+"\t"+"RIF Resistance"+"\t"+"MutCodonNum"+"\t"+"MutCodonVal"+"\t"+"AbsentProbe"+"\t"+"DevelopingProbe"+"\t"+"S1"+"\t"+"S2"+"\t"+"S3"+"\t"+"S4"+"\t"+"S5"+"\t"+"R2"+"\t"+"R4a"+"\t"+"R4b"+"\t"+"R5"+"\n")
  
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
outputList=["1","1","1","1","1","0","0","0","0"]
outputKeys=["S1","S2","S3","S4","S5","R2","R4a","R4b","R5"]

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
					if codon in Nipromut.keys(): #check if the codon is in the mutation list
						Nipromutcodval = sampleCodons[codon]
						ind = sublist.index(mutpos)
						#modify the WT codon to be the new mutated codon
						Nipromutcodval[ind] = altbase
						sampleCodons[codon]=Nipromutcodval					
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
				codnucres=Nipromut[codonPos]
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
			NiproOutputfinal.write(sampleName+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(wtprobeList)+"\t"+",".join(mutprobeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			NiproOutputfinal.write(sampleName+"\t"+"NOT DETECTED\n")
NiproOutputfinal.close()
sys.exit()