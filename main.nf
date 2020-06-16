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

Channel.fromPath(params.rawdata_dir)
  .map { file -> tuple(file.baseName, file) }
  .set { genomes }

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
process get_single_copy {
  label 'busco'
  beforeScript "${params.busco_env}"

  publishDir "${params.outdir}/${params.assembly_completness_dirname}" , mode: 'copy', pattern : "${genome_name}/short_summary*"
  publishDir "${params.outdir}/${params.assembly_completness_dirname}" , mode: 'copy', pattern : "${genome_name}/run_*/*.tsv" , saveAs : { full_table -> "${genome_name}/full_table.tsv" }
  publishDir "${params.outdir}/${params.assembly_completness_dirname}" , mode: 'copy', pattern : "*.faa"

  input:
    set genome_name, file(fasta) from genomes

  output:
    file "${genome_name}/short_summary*" into busco_short_summary
    file "${genome_name}/run_*/full_table.tsv" into busco_full_summary
    file "${genome_name}.sg.faa" into busco_single_copy_proteins

  shell:
    """
    busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${genome_name} -l ${params.odb_path}/${params.odb_name} >& busco.log 2>&1
    catSingleCopyBySpecie.py -f ${genome_name}/run_${params.odb_name}/busco_sequences/single_copy_busco_sequences/ -s ${genome_name} >& catSingleCopyBySpecie.log 2>&1
    """

}

// Concat all single copy genes of all specie into a single file
process concat_single_copy {

  publishDir "${params.outdir}/${params.concatenate_sg_dirname}", mode: 'copy'

  input:
    file '*.faa' from busco_single_copy_proteins.collect()

  output:
    file "single_copy_busco_sequences.fasta" into busco_single_copy_proteins_concat

  shell:
    """
    cat *.faa > single_copy_busco_sequences.fasta
    """
}

// Filter in order to keep at least a minimum of X species which shared a single copy gene (give a list of sequence ID)
process filter_single_copy {
  beforeScript "${params.biopython_env}"

  publishDir "${params.outdir}/${params.shared_sg_dirname}", mode: 'copy'

  input:
    file fasta from busco_single_copy_proteins_concat

  output:
    file "*.faa" into busco_single_copy_proteins_shared

  shell:
    """
    extractSharedSingleCopy.py -f ${fasta} -s ${params.min_species}
    """
}

// Align each orthogroup
process mafft {
  beforeScript "${params.mafft_env}"

  publishDir "${params.outdir}/${params.alignment_dirname}", mode: 'copy'

  input:
    file faa from busco_single_copy_proteins_shared.flatten()

  output:
    file "*.mafft" into busco_single_copy_proteins_toGblocks

  shell:
    """
    mafft --auto ${faa} > ${faa}.mafft
    """
}

// Clean the alignment
process gblocks {
  // As Gblocks exit status is always 1...
  validExitStatus 1
  beforeScript "${params.gblocks_env}"

  publishDir "${params.outdir}/${params.alignment_blocks_dirname}", mode: 'copy'

  input:
    file aln from busco_single_copy_proteins_toGblocks

  output:
    file "*-gb" into busco_single_copy_proteins_toCleanHeader

  script:
    """
    Gblocks ${aln} -t=p -p=n -b3=8 -b4=10 -b5=h
    """
}

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
