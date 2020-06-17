# Introduction

PALLOR (Phylogeny from universAL singLe cOpy oRthologs) is a bioinformatics workflow used to create in a fast way phylogenetic tree from (un)annotated genomes. It can be useful for genome assembly projects when you have your draft assembly and before did any structural annotation process.


Steps of the workflow:
- Make BUSCO analysis on each genome and extract single copy genes (amino acid)
- Concatenate all single copy genes files into a uniq file
- Extract single copy genes shared by a minimal number of species
- Make the alignment (Mafft)
- Clean it (Gblocks)
- Concatenate alignments into a matrix (ElConcatenero)
- Make the tree (IQ-TREE 2)

## Quick start

1. Clone the current gitlab repertory in your working directory on DATARMOR

```
git clone https://gitlab.ifremer.fr/bioinfo/pallor.git
```

2. Once on DATARMOR, edit the conf/params.config file with your analysis parameters or directly using arguments

3. Add a directory with your data (fasta files of genomes)

4. Run the analysis

```
bash pallor.sh
```

## License and Credits
PALLOR is released under the GNU Affero General Public License, Version 3.0. AGPL

It is developped by Alexandre Cormier, bioinformatics engineer at the bioinformatics service of IFREMER (SeBiMER).

-- (c) 2020 - SeBiMER, Ifremer
