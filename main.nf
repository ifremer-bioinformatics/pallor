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
    --min_species [int]                     Keep orthologs presents in a least this minimal number of species (min:2; max: total number of species)

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

  publishDir "${params.outdir}/${params.completness_dirname}" , mode: 'copy', pattern : "${genome_name}/short_summary*"
  publishDir "${params.outdir}/${params.completness_dirname}" , mode: 'copy', pattern : "${genome_name}/run_*/*.tsv" , saveAs : { full_table -> "${genome_name}/full_table.tsv" }
  publishDir "${params.outdir}/${params.completness_dirname}" , mode: 'copy', pattern : "*.faa"

  input:
    set genome_name, file(fasta) from genomes

  output:
    file "${genome_name}/short_summary*" //into completness_short_summary
    file "${genome_name}/run_*/full_table.tsv" //into completness_full_summary
    file "${genome_name}.sg.faa" into single_copy_proteins

  shell:
    """
    busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${genome_name} -l ${params.odb_path}/${params.odb_name} >& busco.log 2>&1
    catSingleCopyBySpecie.py -f ${genome_name}/run_${params.odb_name}/busco_sequences/single_copy_busco_sequences/ -s ${genome_name} >& catSingleCopyBySpecie.log 2>&1
    """

}

// Concat all single copy genes of all specie into a single file
process concat_single_copy {

  publishDir "${params.outdir}/${params.concatenate_dirname}", mode: 'copy'

  input:
    file '*.faa' from single_copy_proteins.collect()

  output:
    file "single_copy_busco_sequences.fasta" into cat_all_single_copy_proteins

  shell:
    """
    cat *.faa > single_copy_busco_sequences.fasta
    """
}

// Filter in order to keep at least a minimum of X species which shared a single copy gene (give a list of sequence ID)
process filter_single_copy {
  beforeScript "${params.biopython_env}"

  publishDir "${params.outdir}/${params.extract_shared_sg_dirname}", mode: 'copy'

  input:
    file fasta from cat_all_single_copy_proteins

  output:
    file "*.faa" into shared_single_copy_proteins

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
    file faa from shared_single_copy_proteins.flatten()

  output:
    file "*.mafft" into aligned_shared_single_copy_proteins

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

  publishDir "${params.outdir}/${params.cleaning_dirname}", mode: 'copy'

  input:
    file aln from aligned_shared_single_copy_proteins

  output:
    file "*-gb" into cleaned_aligned_shared_single_copy_proteins

  script:
    """
    Gblocks ${aln} -t=p -p=n -b3=8 -b4=10 -b5=h
    """
}

// Create the matrix
process concatenation {
  beforeScript "${params.ElConcatenero_env}"

  publishDir "${params.outdir}/${params.matrix_dirname}", mode: 'copy'

  input:
    file '*-gb' from cleaned_aligned_shared_single_copy_proteins.collect()

  output:
    file "*.fas" into concatenated_alignments

  script:
    """
    ElConcatenero.py -if fasta -of fasta -in *-gb -o concatenated
    """
}

process iqtree {
  label 'iqtree'
  beforeScript "${params.iqtree_env}"

  publishDir "${params.outdir}/${params.tree_dirname}", mode: 'copy'

  input:
    file matrix from concatenated_alignments

  output:
    file "*.treefile" into iqtree_tree
    file "*.iqtree" into iqtree_logs

  shell:
    """
    iqtree -s ${matrix} -bb 1000 -alrt 1000 -nt ${task.cpus}
    """
}

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
