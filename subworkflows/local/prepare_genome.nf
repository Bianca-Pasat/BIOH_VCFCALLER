// from nf-core/rnaseq
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BWAMEM_INDEX       } from '../../modules/nf-core/untar/main'
include { BWAMEM2_INDEX                      } from '../../modules/nf-core/bwamem2/index/main'
include { CAT_ADDITIONAL_FASTA              } from '../../modules/local/cat_additional_fasta'
include { SAMTOOLS_FAIDX1 } from '../../modules/nf-core/samtools/faidx'


workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list   : tools to prepare indices for


    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }


    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
   /* ch_gtf = Channel.empty()
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], params.additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, biotype )
        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }
*/


    //
    // Uncompress BWAMEM2 index or generate from scratch if required
    //

    ch_bwamem_index = Channel.empty()
    if ('bwamem2' in prepare_tool_indices) {
        if (params.bwamem_index_provided) {
            if (params.bwamem_index.endsWith('.tar.gz')) {
                ch_bwamem_index = UNTAR_BWAMEM2_INDEX ( [ [:], params.bwamem2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_BWAMEM2_INDEX.out.versions)
            } else {
                ch_bwamem_index = file(params.bwamem2_index)
            }
        } else {
            fasta_meta = (ch_fasta).map{ it -> [[id:it[0].baseName], it] }
            BWAMEM2_INDEX ( fasta_meta )
            ch_bwamem_index = BWAMEM2_INDEX.out.index
            ch_versions     = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }
    }
    
    ch_fai = Channel.empty()
    SAMTOOLS_FAIDX1(ch_fasta)
    ch_fai=SAMTOOLS_FAIDX1.out.fai



    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    meta_fasta            = fasta_meta            //    path: genome.fasta
   // gtf              = ch_gtf              //    path: genome.gtf
    fai              = ch_fai              //    path: genome.fai
    bwamem2_index     = ch_bwamem_index     //    path: bwamem2/index/


    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}
