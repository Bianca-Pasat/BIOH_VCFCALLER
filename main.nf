#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/qualitycontrol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/qualitycontrol

    Website: https://nf-co.re/qualitycontrol
    Slack  : https://nfcore.slack.com/channels/qualitycontrol
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf           = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff')
params.bbsplit_index = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.hisat2_index  = WorkflowMain.getGenomeAttribute(params, 'hisat2')
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { QUALITYCONTROL } from './workflows/qualitycontrol'
include { VCFCALLER } from './workflows/vcfcaller.nf'
include { GLNEXUS } from './modules/local/glnexus/main/'
//include { SRA } from './workflows/sra.nf'
//include { RNAVAR } from './workflows/rnavar.nf'
//include { DENET } from './workflows/denet'
//include { ANNOTATEVCF } from './workflows/annotatevcf.nf'
//include { PREPARE_GENOME } from './subworkflows/local/prepare_genome.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if --input file is empty
ch_input = file(params.input)
//ch_input = file(params.input, checkIfExists: true)
//if (ch_input.isEmpty()) {exit 1, "File provided with --input is empty: ${ch_input.getName()}!"}

// Read in ids from --input file
/*Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:false, sep:'', strip:true)
    .map { it[0] }
    .unique()
    .set { ch_ids }
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Auto-detect input id type
/* def input_type = ''
if (WorkflowMain.isSraId(ch_input, log)) {
    input_type = 'sra'
} else if (WorkflowMain.isSynapseId(ch_input, log)) {
    input_type = 'synapse'
} else {
    exit 1, 'Ids provided via --input not recognised please make sure they are either SRA / ENA / DDBJ or Synapse ids!'
}

if (params.input_type == input_type) {
    if (params.input_type == 'sra') {
        include { SRA } from './workflows/sra'
    } else if (params.input_type == 'synapse') {
        include { SYNAPSE } from './workflows/synapse'
    }
} else {
    exit 1, "Ids auto-detected as ${input_type}. Please provide '--input_type ${input_type}' as a parameter to the pipeline!"
} */
//
// WORKFLOW: Run main nf-core/qualitycontrol analysis pipeline
//
workflow NFCORE_VCF {
VCFCALLER(ch_input)
//GLNEXUS(params.glnexusdir,params.bed,params.glnexusconfig, VCFCALLER.out.sanity)

/*
    if (!params.sra) {
        QUALITYCONTROL (params.input, params.ribo_database_manifest, params.bbsplit_fasta_list)
        ALIGNASSEMBLY (QUALITYCONTROL.out.filtered_reads)


    }else if (params.sra){
    //    ch_ids.collect().view()
        SRA ( ch_ids )
        QUALITYCONTROL(SRA.out.samplesheet,params.ribo_database_manifest, params.bbsplit_fasta_list)
        ALIGNASSEMBLY (QUALITYCONTROL.out.filtered_reads)


}
*/
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_VCF ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
