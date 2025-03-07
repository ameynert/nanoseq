/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it, params.genomes) }
        .map { it -> [ it[0], it[2], it[3], it[4], it[5], it[6], it[1] , it[7] ] }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str ]
def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {
    def meta = [:]
    meta.id  = sample.sample

    // Resolve fasta and gtf file if using iGenomes
    def fasta = false
    def gtf   = false
    if (sample.genome) {
        if (genomeMap && genomeMap.containsKey(sample.genome)) {
            fasta = file(genomeMap[sample.genome].fasta, checkIfExists: true)
            gtf   = file(genomeMap[sample.genome].gtf, checkIfExists: true)
        } else {
            fasta = file(sample.genome, checkIfExists: true)
        }
    }

    // Check if input file and gtf file exists
    input_file = sample.input_file ? file(sample.input_file, checkIfExists: true) : null
    gtf        = sample.gtf        ? file(sample.gtf, checkIfExists: true)        : gtf

    return [ meta, input_file, sample.barcode, fasta, gtf, sample.is_transcripts.toBoolean(), fasta.toString()+';'+gtf.toString() ]
}
