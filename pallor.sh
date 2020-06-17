#!/usr/bin/env bash

. /appli/bioinfo/nextflow/20.01.0/env.sh

#nextflow temp directory
export NXF_TEMP=$SCRATCH/.nxf_temp
mkdir -p $NXF_TEMP

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
#data=/home1/datawork/acormier/sgmm/2-vibrio/0-data/*_R{1,2}*.fastq.gz
#out=/home1/datawork/acormier/sgmm/2-vibrio/celia/

nextflow -DXms=2G -DXmx=8G -trace nextflow.executor run main.nf --min_species 4 

#deactivate nextflow environment
. /appli/bioinfo/nextflow/20.01.0/delenv.sh

###
date
###
