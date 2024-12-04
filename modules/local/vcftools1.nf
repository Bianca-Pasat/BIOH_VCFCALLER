process VCFTOOLS1{
    
    tag "${base}"
    
    //publishDir "$params.outdir/vcftools1/", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/vcf:vcftools'}"
    input:
    
    path(vcf_vcftools) //from ch_out_vcf
    val(pass) //from params.pass
    val(recode) // from params.recode
    val(recode_INFO_all) // from params.recode_INFO_all
    val(minall) // from params.minall
    val(maxall) // from params.maxall
    val(gq) // from params.gq
    val(dp) // from params.dp
     
    output:
    path("*.vcf"), emit: ch_out_vcftools
//    path("*_vcftools_copy2.vcf"),             emit: ch_out_vcftools_copy
    
    script:
    """
    vcftools --vcf $vcf_vcftools $pass $recode $recode_INFO_all $minall $maxall $gq $dp --out ${vcf_vcftools.baseName}_recoded.vcf 
    """
}

//&& cp *_recoded.vcf ${vcf_vcftools.baseName}_vcftools_copy2.vcf
