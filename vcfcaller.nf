//nextflow.enable.dsl=2

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowDenet.initialise(params, log)
/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { ALIGNASSEMBLY } from './alignassembly.nf'
include { SORTBAM } from '../subworkflows/local/sortbam.nf'
include { INPUT_CHECK } from '../subworkflows/local/input_check.nf'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { GUNZIP          } from '../modules/nf-core/gunzip/main'

include { FASTP                  } from '../modules/nf-core/fastp/main'

include { DEEPVARIANT1   } from '../modules/nf-core/deepvariant1/main'

ch_versions = Channel.empty() 
workflow VCFCALLER{

    take:
    ch_input

    main:
    INPUT_CHECK (
        ch_input
    ).reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))
    
    FASTP(ch_cat_fastq,  params.fastp_qualified_quality_phred, params.fastp_unqualified_percent_limit,params.fastp_trim_poly_g, params.fastp_length_required,false, false,false)

    // QUALITYCONTROL(ch_cat_fastq)
    
    //GUNZIP(FASTP.out.reads)
    //ch_unzipped = GUNZIP.out.gunzip
    
    //ch_unziped = QUALIYTYCONTROL.out.filtered_reads)
    
    ALIGNASSEMBLY(FASTP.out.reads)
    /// ATTENTION I CHANGED TO SAM AND SAMTOOLS COLLATE HAS ALSO SAM!!!!!!

    /// TO DO ATTENTION I NEED TO CONVERT BACK TO TUPLE AAAAAALLLL THE SAMTOOLS PROCESSES!!!!!!!!!!!!!!!!!!!!!!!1
    SORTBAM(ALIGNASSEMBLY.out.genome_sam)

   /* bam_channel = Channel.fromPath(params.example_bam)
    bam_tuple_channel = bam_channel.map { bam_file -> tuple(bam_file) } // NOT WORKING!!!!!!!!!
    SORTBAM(bam_tuple_channel)
*/
    //DEEPVARIANT1(SORTBAM.out.bam_only,ALIGNASSEMBLY.out.fasta, ALIGNASSEMBLY.out.fai)
    DEEPVARIANT1(params.testbam,params.testbai,params.fasta1, params.fai1)
    //GLNEXUS(params.deepvariant_outdir, DEEPVARIANT.out.sanitycheck,params.glnexus_config)
    emit:
    sanity = DEEPVARIANT1.out.versions
  }

