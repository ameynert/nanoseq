import org.yaml.snakeyaml.Yaml

/*
 * This file holds several functions used to perform standard checks for the nf-core pipeline template.
 */

class Checks {

    static void check_conda_channels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        def required_channels = ['conda-forge', 'bioconda', 'defaults']
        def conda_check_failed = !required_channels.every { ch -> ch in channels }

        // Check that they are in the right order
        conda_check_failed |= !(channels.indexOf('conda-forge') < channels.indexOf('bioconda'))
        conda_check_failed |= !(channels.indexOf('bioconda') < channels.indexOf('defaults'))

        if (conda_check_failed) {
            log.warn "=============================================================================\n" +
            "  There is a problem with your Conda configuration!\n\n" +
            "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
            "  Please refer to https://bioconda.github.io/user/install.html#set-up-channels\n" +
            "  NB: The order of the channels matters!\n" +
            "==================================================================================="
        }
    }

    static void aws_batch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            assert (params.awsqueue && params.awsregion) : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            // related: https://github.com/nextflow-io/nextflow/issues/813
            assert params.outdir.startsWith('s3:')       : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
            // Prevent trace files to be stored on S3 since S3 does not support rolling files.
            assert !params.tracedir.startsWith('s3:')    :  "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
        }
    }

    static void hostname(workflow, params, log) {
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (params.hostnames) {
            def hostname = "hostname".execute().text.trim()
            params.hostnames.each { prof, hnames ->
                hnames.each { hname ->
                    if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                        log.info "=${colors.yellow}====================================================${colors.reset}=\n" +
                        "${colors.yellow}WARN: You are running with `-profile $workflow.profile`\n" +
                        "      but your machine hostname is ${colors.white}'$hostname'${colors.reset}.\n" +
                        "      ${colors.yellow_bold}Please use `-profile $prof${colors.reset}`\n" +
                        "=${colors.yellow}====================================================${colors.reset}="
                    }
                }
            }
        }
    }

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The pipeline\n" +
        "  https://doi.org/10.5281/zenodo.3697959\n\n" +
        "* The nf-core framework\n" +
        "  https://dx.doi.org/10.1038/s41587-020-0439-x\n" +
        "  https://rdcu.be/b1GjZ\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    // Print a warning if using GRCh38 assembly from igenomes.config
    static void ncbi_genome_warn(log) {
        log.warn "=============================================================================\n" +
        "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
        "  Auto-activating '--skip_biotype_qc' parameter to circumvent the issue below:\n" +
        "  https://github.com/nf-core/rnaseq/issues/460.\n\n" +
        "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
        "==================================================================================="
    }

    // Print a warning if using a UCSC assembly from igenomes.config
    static void ucsc_genome_warn(log) {
        log.warn "=============================================================================\n" +
        "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
        "  Auto-activating '--skip_biotype_qc' parameter to circumvent the issue below:\n" +
        "  https://github.com/nf-core/rnaseq/issues/460.\n\n" +
        "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
        "==================================================================================="
    }

    // Print a warning if --skip_alignment has been provided
    static void skip_alignment_warn(log) {
        log.warn "=============================================================================\n" +
        "  '--skip_alignment' parameter has been provided.\n" +
        "  Skipping alignment, quantification and all downstream QC processes.\n" +
        "==================================================================================="
    }

}
