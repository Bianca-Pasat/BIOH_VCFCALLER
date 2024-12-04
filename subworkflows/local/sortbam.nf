
// runs SAMPLE_QC from either reads or bam files
// or both after alignment

include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP   } from '../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_COLLATE    } from '../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FIXMATE  } from '../../modules/nf-core/samtools/fixmate/main'

/// TO DO ATTENTION I NEED TO CONVERT BACK TO TUPLE AAAAAALLLL THE SAMTOOLS PROCESSES!!!!!!!!!!!!!!!!!!!!!!!1
workflow SORTBAM {

    take:
    bam        // channel: [mandatory] [ val(meta), path(bam) ]

    main:
    ch_versions = Channel.empty()

    // samtools stats block needs the bam file to be sorted
    // and indexed

    SAMTOOLS_COLLATE ( bam )
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions.first())

    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions.first())

    // additionally, the modules require a single channel containing
    // both the bam file and its index
    // which is the reason for creating an additional channel that joins
    // both of them:

    SAMTOOLS_SORT(SAMTOOLS_FIXMATE.out.bam)

    SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam) 

    emit:
    bam_only = SAMTOOLS_MARKDUP.out.bam  // channel: [ val(meta), path(bam) ]
    //bam_bai  = bam_bai                // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions            // channel: [ versions.yml ]
}


