// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0"
    } else {
        container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"
    }

    input:
    tuple val(meta), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)

    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam
    path "versions.yml"           , emit: versions

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
    def junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf $gtf" : ""
    """
    graphmap2 \\
        align \\
        $preset \\
        $junctions \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $fastq \\
        -o ${meta.id}.sam \\
        --extcigar

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
