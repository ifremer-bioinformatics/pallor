#!/usr/bin/env python

import argparse

def getArgs():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-f',dest="fasta",type=argparse.FileType('r'),required=True,help='')
	
	arg = parser.parse_args()
	
	return arg

def main(args):
	
	for line in args.fasta:
		if line.startswith('>'):
			eog, specie, contig, coordinates = line.split(':')
			print('>'+'.'.join(specie.split('.')[:-1]))
		else:
			print(line[:-1])

if __name__ == '__main__':
	args = getArgs()
	main(args)
