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

## Dependencies

- Nextflow >= 20.04.x
- Singularity or Docker

## Quick start

1. Clone the current gitlab repertory 

```
git clone https://gitlab.ifremer.fr/bioinfo/pallor.git
```

2. Add a directory with your data (fasta files of genomes with ```.fna``` extension). It should be look like:

3. Run the analysis

```
nextflow run main.nf -profile custom,singularity
```

## Parameters

- ```--rawdata_dir```: Path to input directory with raw data files with ```.fna``` extension [path]
- ```--name```: Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic [str]
- ```--projectName```: Name of the project being analyzed [str]
- ```--odb_path```: Path to all BUSCO ODB databases [path]
- ```--odb_name```: Specify the name of the BUSCO lineage to be used [str]
- ```--min_species```: Keep orthologs presents in a least ```--min_species``` minimal number of species (min:2; max: total number of species) [int]

## Examples

- 1 - Run using Singularity on 20 microsporidia species

```
nextflow run main.nf -profile base,singularity --projectName run-test --rawdata_dir test_data/ --min_species 18 --odb_path /home/ref-bioinfo/tools_data/busco/v4 --odb_name microsporidia_odb10
```

- 2 - Run using Docker and parameters set in conf/custom.conf

```
nextflow run main.nf -profile custom,docker
```

- 3 - Same with HPC configuration

```
nextflow run main.nf -profile custom,docker -c /appli/bioinfo/hpc/nextflow/ifremer.config
```

## License and Credits
PALLOR is released under the GNU Affero General Public License, Version 3.0. AGPL

It is developped by Alexandre Cormier, bioinformatics engineer at the bioinformatics service of IFREMER (SeBiMER).

-- (c) 2020 - SeBiMER, Ifremer
