#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

def getArgs():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-f',dest="fasta", type=str,required=True,help='Fasta file with all single copy busco sequences')
	parser.add_argument('-s',dest="species",type=int,required=True,help='Minimal number of species per cluster')
	
	args = parser.parse_args()
	
	return args

def main(args):
	
	seqBUSCO = {}
	
	for fasta in SeqIO.parse(args.fasta, "fasta"):
		buscoID, specie, contig, coordinates = fasta.id.split(':')
		if not buscoID in seqBUSCO:
			seqBUSCO[buscoID] = {'seq':[fasta.id], 'species':[specie]}
		else:
			seqBUSCO[buscoID]['seq'].append(fasta.id)
			seqBUSCO[buscoID]['species'].append(specie)
	
	
	for seqID, infos in sorted(seqBUSCO.items()):
		if (len(infos['species']) - len(set(infos['species']))) == 0 and len(infos['species']) >= args.species:
			outpout = open(seqID+'.lst','w')
			for seq in infos['seq']:
				outpout.write(seq+'\n')
			outpout.close()
		
if __name__ == '__main__':
	args = getArgs()
	main(args)