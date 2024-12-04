process GLNEXUS {
 //   tag "$txt"
    //tag "$meta.id"
    //tag "$names"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/glnexus:latest'}"



    input:


    //tuple val(meta), path(txt, stageAs: "*")
    //tuple val(meta), path(txt, stageAs: "input*/*")
    path(pathi)
   
    path(bed)
    path(config)
    path(sanity)
    //val(names)



    output:
    path "*.bcf"        , emit: mergedBcf
//    path "test"        , emit: test
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when
    def vcfs = "find -L $pathi -name '*.g.vcf.gz'"
    /*if (!vcfs) {
        error "No VCF files found in directory: ${pathi}"
    }
*/
    script:
    

    
    """
    find -L $pathi -name '*.g.vcf.gz' | xargs /glnexus_cli --config $config --bed $bed > merged.bcf






    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus
    END_VERSIONS
    """
}
