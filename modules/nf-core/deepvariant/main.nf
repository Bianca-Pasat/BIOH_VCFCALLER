process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    // FIXME Conda is not supported at the moment
    // BUG https://github.com/nf-core/modules/issues/1754
    // BUG https://github.com/bioconda/bioconda-recipes/issues/30310
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'google/deepvariant:latest'}"

    input:
    tuple val(meta), path(input) // path(index) //path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    //tuple val(meta4), path(gzi)
    //tuple val(meta5), path(par_bed)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi
    path "wooowheeee" , emit: sanitycheck
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions=${intervals}" : ""
    //def par_regions = par_bed ? "--par_regions_bed=${par_bed}" : ""
    // WARN https://github.com/nf-core/modules/pull/5801#issuecomment-2194293755
    // FIXME Revert this on next version bump
    def VERSION = '1.6.1'

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --include_med_dp
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        ${args} \\
        --intermediate_results_dir=tmp \\
        --model_type=WGS" \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    touch wooowheeee
    
    "${task.process}":
        deepvariant: $VERSION
    END_VERSIONS
    """
}