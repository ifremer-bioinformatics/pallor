#!/usr/bin/env nextflow

/*
========================================================================================
                          PALLOR: Phylogeny from universAL singLe cOpy oRthologs
========================================================================================
 CELIA Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.ifremer.fr/bioinfo/pallor
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  // Add to this help message with new command line parameters
  log.info SeBiMERHeader()
  log.info"""
  Usage:

  The typical command for running the pipeline after filling the conf/base.config file is as follows :
    nextflow run main.nf

    Mandatory arguments:
    --rawdata_dir [path]                    Path to input directory with raw data files

    Other options:
    --outdir [path]                         The output directory where the results will be saved
    -w/--work-dir                           The temporary directory where intermediate data will be saved
    -name [str]                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --projectName [str]                     Name of the project being analyzed

    Single copy orthologs:
    --odb_path [path]                       Path to all BUSCO ODB databases
    --odb_name [str]                        Specify the name of the BUSCO lineage to be used

  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

/*
* SET UP CONFIGURATION VARIABLES
*/

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

Channel.fromPath(params.indir)
  .map { file -> tuple(file.baseName, file) }
  .into { genomes }

/*
* PIPELINE INFO
*/
// Header log info
log.info SeBiMERHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Project Name']     = params.projectName
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName

// if (params.email || params.email_on_fail) {
//   summary['E-mail Address']    = params.email
//   summary['E-mail on failure'] = params.email_on_fail
// }

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Launch BUSCO and cat all single copy genes into a single file for each specie
process busco {
  label 'busco'
  beforeScript "${params.busco_env}"

  publishDir "${params.outdir}/${params.busco_dirname}" , mode: 'copy'

  input:
    set name, file(fasta) from genomes

  output:
    file ("run_${name}/*.{txt,tsv}") into busco_summary_results
    file ("run_${name}/single_copy_busco_sequences_${name}.faa") into busco_single_copy_proteins

  shell:
    """
    run_BUSCO.py -i ${fasta} -o ${name} -m genome -l ${BUSCOodb} --cpu ${params.busco.cpus}
    cat run_${name}/single_copy_busco_sequences/*.faa > run_${name}/single_copy_busco_sequences_${name}.faa
    """
}

// Concat all single copy genes of all specie into a single file
// process concat_busco {
//
//   publishDir "${params.outdir}/1-busco", mode: 'copy'
//
//   input:
//     file fasta from busco_single_copy_proteins.collect()
//
//   output:
//     file "single_copy_busco_sequences.fasta" into busco_single_copy_proteins_concat
//
//   shell:
//     """
//     cat ${fasta} > single_copy_busco_sequences.fasta
//     """
// }

// Filter in order to keep at least a minimum of X species which shared a single copy gene (give a list of sequence ID)
// process filter_single_copy {
//   beforeScript "${params.biopython_env}"
//
//   publishDir "${params.outdir}/2-single-copy", mode: 'copy'
//
//   input:
//     file fasta from busco_single_copy_proteins_concat
//
//   output:
//     file "*.lst" into busco_single_copy_proteins_list //mode flatten
//     // flatten to allow parallelization at the next step
//
//   // add a limitation on -s option in order than it can be superior to the number of input fasta
//   shell:
//     """
//     single_copies_busco.py -f ${fasta} -s params.species
//     """
// }

// Extract each ortholgous sequences
// process seqtk {
//   beforeScript "${params.seqtk_env}"
//
//   publishDir "${params.outdir}/2-single-copy", mode: 'copy'
//
//   input:
//     file lst from busco_single_copy_proteins_list.flatten()
//     file fasta from busco_single_copy_proteins_concat
//
//   output:
//     file "*.faa" into busco_single_copy_proteins_toMafft
//
//   shell:
//     """
//     seqtk subseq -l 70 ${fasta} ${lst} > ${lst}.faa
//     """
// }

// Align each orthogroup
// process mafft {
//   beforeScript "${params.mafft_env}"
//
//   publishDir "${params.outdir}/3-mafft", mode: 'copy'
//
//   input:
//     file faa from busco_single_copy_proteins_toMafft
//
//   output:
//     file "*.mafft" into busco_single_copy_proteins_toGblocks
//
//   shell:
//     """
//     mafft --auto ${faa} > ${faa}.mafft
//     """
// }

// Clean the alignment
// process gblocks {
//   // As Gblocks exit status is always 1...
//   validExitStatus 1
//   beforeScript "${params.gblocks_env}"
//
//   publishDir "${params.outdir}/4-gblocks", mode: 'copy'
//
//   input:
//     file aln from busco_single_copy_proteins_toGblocks
//
//   output:
//     file "*-gb" into busco_single_copy_proteins_toCleanHeader
//
//   script:
//     """
//     Gblocks ${aln} -t=p -p=n -b3=8 -b4=10 -b5=h
//     """
// }

// Clean the header for higher readability in the final tree
// Default: EOG093001PG:NGRA.fna:NGRA000074:9808-10725
// Cleaned: NGRA
// process cleanHeader {
//   beforeScript "${params.biopython_env}"
//
//   publishDir "${params.outdir}/5-cat", mode: 'copy'
//
//   input:
//     file gblocks from busco_single_copy_proteins_toCleanHeader
//
//   output:
//     file "*.fas" into busco_single_copy_proteins_toFASconCAT
//
//   script:
//     """
//     cleanHeaderFromBuscoSG.py -f ${gblocks} > ${gblocks}.fas
//     """
// }

// Create the matrix
// process FASconCAT {
//   beforeScript "${params.fasconcat_env}"
//
//   publishDir "${params.outdir}/5-cat", mode: 'copy'
//
//   input:
//     file gblocks from busco_single_copy_proteins_toFASconCAT.collect()
//
//   output:
//     file "FcC_smatrix.fas" into busco_single_copy_proteins_toProTest
//
//   script:
//     """
//     FASconCAT_v1.0.pl -s
//     """
// }

// Find the better paramters for RAxML
// process iqtree {
//   beforeScript "${params.iq-tree_env}"
//
//   publishDir "${params.outdir}/5-tree", mode: 'copy'
//
//   input:
//     file matrix from busco_single_copy_proteins_toProTest
//
//   output:
//     file "*.newick" into final_tree
//
//   script:
//     """
//     iqtree -s ${matrix} -bb 1000 -alrt 1000 -nt ${task.cpus}
//     """
// }

/* Other functions */
def SeBiMERHeader() {
    // Log colors ANSI codes
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_cyan}--------------------------------------------------${c_reset}-
    ${c_blue}    __  __  __  .       __  __  ${c_reset}
    ${c_blue}   \\   |_  |__) | |\\/| |_  |__)  ${c_reset}
    ${c_blue}  __\\  |__ |__) | |  | |__ |  \\  ${c_reset}
                                            ${c_reset}
    ${c_yellow}  PALLOR: Phylogeny from universAL singLe cOpy oRthologs${c_reset}
    -${c_cyan}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
