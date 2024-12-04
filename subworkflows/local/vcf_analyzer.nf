//
// ANNOTATION
//

include { SNPEFF_ANNOTATE } from '../../modules/local/snpeffann'
include { SNPSIFT_ANNOTATE } from '../../modules/local/snpsift_annotate'
include { SNPSIFT_EXTRACT_INDELS } from '../../modules/local/snpsift_extract_indels'
include { SNPSIFT_EXTRACT_SNPS } from '../../modules/local/snpsift_extract_snps'
include { SNPSIFTFILTER as SNPSIFT_FILTER } from '../../modules/local/snpsiftfilter'
include { SNPSIFT_FILTER as SNPSIFT_FILTER2 } from '../../modules/local/snpsift_filter'
include { SNPSIFT_FULLY_ANNOTATE } from '../../modules/local/snpsift_fully_annotate'
include { SPLIT_VARIANTS1 } from '../../modules/local/split_variants1'
include { SPLIT_VARIANTS2 } from '../../modules/local/split_variants2'
include { VCFTOOLS1 } from '../../modules/local/vcftools1'
//include { ENSEMBLVEP_ANNOTATE } from '../../subworkflows/nf-core/ensemblvep_annotate'
//include { INPUTVALIDATOR } from '../../modules/local/inputvalidator'

workflow VCF_ANALYZER {
    take:
    vcf          // channel: [ val(meta), vcf ]
    //tools
    //snpeff_db
    //snpeff_cache



    /*vep_genome
    vep_species
    vep_cache_version
    vep_cache
*/
    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    //ch_vcf_ann  = Channel.empty()


    VCFTOOLS1(vcf, params.vcfpass, params.vcfrecode,params.vcf_recode_info_all, params.vcf_minall, params.vcf_maxall, params.vcf_gq, params.vcf_dp)
    SNPEFF_ANNOTATE(VCFTOOLS1.out.ch_out_vcftools, params.assembly, params.snpeffdb)
    SNPSIFT_FILTER(SNPEFF_ANNOTATE.out.ch_snpeff_out, params.effect_expression)
    //(ch_split_var1,ch_split_var2)=SNPSIFT_FILTER.out.ch_snpsift_filter_out.into(2)
    SPLIT_VARIANTS1(SNPSIFT_FILTER.out.ch_snpsift_filter_out, params.vcfpass, params.vcfrecode,params.vcf_recode_info_all, params.vcf_minall, params.vcf_maxall, params.vcf_gq, params.vcf_dp, params.keep_indels)
    SPLIT_VARIANTS2(SNPSIFT_FILTER.out.ch_snpsift_filter_out, params.vcfpass, params.vcfrecode,params.vcf_recode_info_all, params.vcf_minall, params.vcf_maxall, params.vcf_gq, params.vcf_dp, params.keep_snps)
    SNPSIFT_ANNOTATE(SPLIT_VARIANTS1.out.ch_indels_out,SPLIT_VARIANTS2.out.ch_snps_out, params.snpsiftdb1, params.snpsiftdb2, params.snpsiftdb1_index, params.snpsiftdb2_index)
    SNPSIFT_FULLY_ANNOTATE(SNPSIFT_ANNOTATE.out.ch_snpsift_anno_indels_out, SNPSIFT_ANNOTATE.out.ch_snpsift_anno_snps_out, params.cosmicdb, params.cosmicdb_index)
    //(ch_fulanno_snps,ch_fulanno_snps_bio)=SNPSIFT_FULLY_ANNOTATE.out.ch_snpsift_full_anno_snps_out.into(2)
    //(ch_fulanno_indels,ch_fulanno_indels_bio)=SNPSIFT_FULLY_ANNOTATE.out.ch_snpsift_full_anno_indels_out.into(2)
    SNPSIFT_FILTER2(SNPSIFT_FULLY_ANNOTATE.out.ch_snpsift_full_anno_indels_out, SNPSIFT_FULLY_ANNOTATE.out.ch_snpsift_full_anno_snps_out, params.snpsexpr, params.indelsexpr)
    //(ch_sig_snps_log_in,ch_sig_snps_bio,ch_sig_snps)=SNPSIFT_FILTER2.out.ch_snpsift_sig_snps_out.into(3)
    //(ch_sig_indels_log_in,ch_sig_indels_bio,ch_sig_indels)=SNPSIFT_FILTER2.out.ch_snpsift_sig_indels_out.into(3)
    SNPSIFT_EXTRACT_INDELS(SNPSIFT_FILTER2.out.ch_snpsift_sig_indels_out)
    SNPSIFT_EXTRACT_SNPS(SNPSIFT_FILTER2.out.ch_snpsift_sig_snps_out)
    

/*
    if (tools.contains('snpeff') || tools.contains('merge')) {
        SNPEFF_ANNOTATE (
            vcf,
            snpeff_db,
            snpeff_cache
        )
        ch_vcf_ann  = ch_vcf_ann.mix(SNPEFF_ANNOTATE.out.vcf_tbi)
        ch_reports  = ch_reports.mix(SNPEFF_ANNOTATE.out.reports)
        ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions.first())
    }

    if (tools.contains('merge')) {
        vcf_ann_for_merge = SNPEFF_ANNOTATE.out.vcf_tbi.map{ meta, vcf, tbi -> [meta, vcf] }
        MERGE_ANNOTATE (
            vcf_ann_for_merge,
            vep_genome,
            vep_species,
            vep_cache_version,
            vep_cache
        )
        ch_vcf_ann  = ch_vcf_ann.mix(MERGE_ANNOTATE.out.vcf_tbi)
        ch_reports  = ch_reports.mix(MERGE_ANNOTATE.out.reports)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE.out.versions.first())
    }

    if (tools.contains('vep')) {
        ENSEMBLVEP_ANNOTATE (
            vcf,
            vep_genome,
            vep_species,
            vep_cache_version,
            vep_cache
        )
        ch_vcf_ann  = ch_vcf_ann.mix(ENSEMBLVEP_ANNOTATE.out.vcf_tbi)
        ch_reports  = ch_reports.mix(ENSEMBLVEP_ANNOTATE.out.reports)
        ch_versions = ch_versions.mix(ENSEMBLVEP_ANNOTATE.out.versions.first())
    }
    */

    emit:
        out_indels = SNPSIFT_EXTRACT_INDELS.out.ch_extracted_indels
        out_snps = SNPSIFT_EXTRACT_SNPS.out.ch_extracted_snps
        reports     = ch_reports
        versions    = ch_versions
}
