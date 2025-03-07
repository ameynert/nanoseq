// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

//params.options = [:]
//def options    = initOptions(params.options)

process BEDTOOLS_GENOMECOV {
    label 'process_medium'

    conda     (params.enable_conda ? "bioconda::bedtools=2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    }

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(bam), path(bai)

    output:
    tuple val(meta), path(sizes), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"                             , emit: versions

    script:
    split = (params.protocol == 'DNA' || is_transcripts) ? "" : "-split"
    """
    bedtools \\
        genomecov \\
        -split \\
        -ibam ${bam[0]} \\
        -bg \\
        | bedtools sort > ${meta.id}.bedGraph
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
