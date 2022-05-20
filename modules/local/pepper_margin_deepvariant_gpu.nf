process PEPPER_MARGIN_DEEPVARIANT_GPU {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/kishwars/pepper_deepvariant:r0.8-gpu' :
        'docker.io/kishwars/pepper_deepvariant:r0.8-gpu' }"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*vcf.gz")     ,  emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") ,  emit: tbi
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ""
    prefix      = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions ${intervals}" : ""
    //def gvcf    = params.make_gvcf ? "--gvcf" : ""

    """
    mkdir -p "${prefix}"

    run_pepper_margin_deepvariant call_variant \\
        -b "${input}" \\
        -f "${fasta}" \\
        -o "." \\
        -p "${prefix}" \\
        -t ${task.cpus} \\
        ${regions} \\
        -gpu \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_margin_deepvariant: \$(run_pepper_margin_deepvariant --version | sed 's/VERSION: //g')
    END_VERSIONS
    """
}
