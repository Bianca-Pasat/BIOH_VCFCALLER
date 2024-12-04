/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    aligners       : ['bwamem2'],
]


/*
// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
//ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { BWAMEM2_MEM } from '../modules/nf-core/bwamem2/mem/main/'
//include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as DEDUP_UMI_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
//include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as DEDUP_UMI_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
//include { ALIGN_HISAT2 } from '../subworkflows/nf-core/align_hisat2'
//include { BAM_MARKDUPLICATES_PICARD     } from '../subworkflows/nf-core/bam_markduplicates_picard'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
include { MULTIQC                     } from '../modules/local/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { STRINGTIE_STRINGTIE } from '../modules/nf-core/stringtie/stringtie/main'
//include { RANKPROD } from '../modules/local/rankprod/main'
include { CONCETRATE } from '../modules/local/concetrate/main'
include { CUT } from '../modules/local/cut/main'
//include { IDS_TO_STRING } from '../modules/local/ids_to_string/main'
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/





workflow ALIGNASSEMBLY {


    take:
    ch_input       // channel: [ val(meta), [ reads ] ]

    main:


    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        params.toolindex, 
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    // BWAMEM2 requires unziped fasta
    //BWAMEM2_MEM(ch_input, PREPARE_GENOME.out.bwamem2_index, false) //true: sort bam
    
    // sets bwamem2index to correct input
    //bwamem2index  = params.fasta ? params.bwamem_index  ? Channel.fromPath(params.bwamem_index).collect()  : PREPARE_GENOME.out.bwamem2_index : []

    
    //here  last comment
    BWAMEM2_MEM(ch_input, params.bwamem_index, params.bwamem_indexpath,false) //true: sort bam
    
    
    
    //BWAMEM2_MEM(ch_input, PREPARE_GENOME.out.bwamem2_index, PREPARE_GENOME.out.meta_fasta,true) //true: sort bam
    
    //ch_genome_bam = BWAMEM2_MEM.out.bam
    //ch_genome_sam = BWAMEM2_MEM.out.sam


    ch_genome_bam = Channel.empty()
    ch_genome_sam = Channel.empty()



/*
        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            DEDUP_UMI_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0])
            )
            ch_genome_bam        = DEDUP_UMI_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = DEDUP_UMI_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = DEDUP_UMI_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = DEDUP_UMI_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index = DEDUP_UMI_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS_GENOME.out.versions)
        }




    //
    // MODULE: Run Preseq
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_markduplicates) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam
        )
        ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
        if (params.bam_csi_index) {
            ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }



*/



emit:
    genome_bam = ch_genome_bam
    genome_sam = ch_genome_sam

    fasta = PREPARE_GENOME.out.fasta
    fai = PREPARE_GENOME.out.fai
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

    //
    // MODULE: Pipeline reporting
    //
    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
*/

