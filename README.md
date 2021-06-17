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

2. Add a directory with your data (fasta files of genomes with ```.fna``` extension). It should be look like:
```
pallor/
├── bin
│   ├── catSingleCopyBySpecie.py
│   └── extractSharedSingleCopy.py
├── conf
│   ├── base.config
│   ├── custom.config
│   ├── reports.config
│   └── resources.config
├── containers
│   └── pallor
│       ├── Dockerfile
│       └── environment.yml
├── LICENSE
├── main.nf
├── nextflow.config
├── pallor.sh
├── README.md
└── test_data
    ├── GCF_000091225.1_ASM9122v1_genomic.fna
    ├── GCF_000146465.1_ASM14646v1_genomic.fna
    ├── GCF_000277815.2_ASM27781v3_genomic.fna
    └── GCF_000280035.1_ASM28003v2_genomic.fna
```

3. Once, edit the conf/custom.config file with your analysis parameters or directly using arguments

4. Run the analysis

```
nextflow run main.nf -profile custom,<docker/singularity>
```

## License and Credits
PALLOR is released under the GNU Affero General Public License, Version 3.0. AGPL

It is developped by Alexandre Cormier, bioinformatics engineer at the bioinformatics service of IFREMER (SeBiMER).

-- (c) 2020 - SeBiMER, Ifremer
