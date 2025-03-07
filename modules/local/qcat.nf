// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QCAT {
    tag "$input_path"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:input_path) }

    conda     (params.enable_conda ? "bioconda::qcat=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/qcat:1.1.0--py_0"
    } else {
        container "quay.io/biocontainers/qcat:1.1.0--py_0"
    }

    input:
    path input_path

    output:
    path "fastq/*.fastq.gz"               , emit: fastq
    path "versions.yml"                   , emit: versions

    script:
    def detect_middle = params.qcat_detect_middle ? "--detect-middle $params.qcat_detect_middle" : ""
    """
    ## Unzip fastq file
    ## qcat doesnt support zipped files yet
    FILE=$input_path
    if [[ \$FILE == *.gz ]]
    then
    zcat $input_path > unzipped.fastq
    FILE=unzipped.fastq
    fi
    qcat  \\
    -f \$FILE \\
    -b ./fastq \\
    --kit $params.barcode_kit \\
    --min-score $params.qcat_min_score \\
    $detect_middle

    ## Zip fastq files (cannot find pigz command)
    gzip fastq/*

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(qcat --version 2>&1 | sed 's/^.*qcat //; s/ .*\$//')
    END_VERSIONS
    """
}
