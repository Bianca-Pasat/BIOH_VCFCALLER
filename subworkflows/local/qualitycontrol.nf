


def getTrimGaloreReadsAfterFiltering(log_file) {
        def total_reads = 0
        def filtered_reads = 0
        log_file.eachLine { line ->
            def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
            def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
            if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
            if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
        return total_reads - filtered_reads
    }



include { FASTP } from '../modules/local/fastp/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '/subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


workflow QUALITYCONTROL {

    take:
    input // reads
    main:



    ch_versions = Channel.empty()

    FASTP(ch_cat_fastq, false, false,true, params.fastp_unqualified_percent_limit,params.fastp_unqualified_percent_limit, params.fastp_trim_poly_g, params.fastp_length_required)

   //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_umi_extract,
        params.skip_trimming,
        params.umi_discard_read
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)



    //
    // Filter channels to get samples that passed minimum trimmed read count
    //

    ch_fail_trimming_multiqc = Channel.empty()
    ch_filtered_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    if (!params.skip_trimming) {
        ch_filtered_reads
            .join(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .map {
                meta, reads, trim_log ->
                    if (!meta.single_end) {
                        trim_log = trim_log[-1]
                    }
                    num_reads =  getTrimGaloreReadsAfterFiltering(trim_log)

                    [ meta, reads, num_reads ]
            }
            .set { ch_num_trimmed_reads  }

        ch_num_trimmed_reads
            .map { meta, reads, num_reads -> if (num_reads > params.min_trimmed_reads) [ meta, reads ] }
            .set { ch_filtered_reads }

        ch_num_trimmed_reads
            .map {
                meta, reads, num_reads ->
                if (num_reads <= params.min_trimmed_reads) {
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .set { ch_num_trimmed_reads }


    }





    //


emit:
filtered_reads = ch_filtered_reads
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
