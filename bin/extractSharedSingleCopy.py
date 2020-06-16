#!/usr/bin/env python

import argparse
from Bio import SeqIO

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f',dest="fasta", type=str,required=True,help='Fasta file with all single copy busco sequences')
    parser.add_argument('-s',dest="species",type=int,required=True,help='Minimal number of species per cluster')

    args = parser.parse_args()

    return args

def main(args):

    sg_odb = {}

    # 1 - Store each sequence of each specie by orthologs
    for fasta in SeqIO.parse(args.fasta, "fasta"):
        og_id, specie, contig, coordinates = fasta.id.split(':')
        if not og_id in sg_odb:
            sg_odb[og_id] = { specie: [fasta.id, fasta.seq]}
        else:
            sg_odb[og_id][specie] = [fasta.id, fasta.seq]

    # 2 - Pass on each ortholgous sequences and write if shared by a minimal number of species
    for og_id, species in sg_odb.items():
        if len(species) >= args.species:
            og_file = open(og_id+'.faa','w')
            for specie in species:
                # header = species[specie][0]
                seq = str(species[specie][1])
                og_file.write('>'+specie+'\n')
                while len(seq) > 0:
                    og_file.write(seq[:70]+'\n')
                    seq = seq[70:]

if __name__ == '__main__':
    args = getArgs()
    main(args)
