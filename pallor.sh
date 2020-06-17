#!/usr/bin/env bash

. /appli/bioinfo/nextflow/20.01.0/env.sh

#nextflow temp directory
export NXF_TEMP=$SCRATCH/.nxf_temp
mkdir -p $NXF_TEMP

# Launch pipeline
nextflow -DXms=2G -DXmx=8G -trace nextflow.executor run main.nf --min_species 4

#deactivate nextflow environment
. /appli/bioinfo/nextflow/20.01.0/delenv.sh

###
date
###
