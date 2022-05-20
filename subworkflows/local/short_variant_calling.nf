/*
 * Short variant calling test
 */

include { MEDAKA_VARIANT                        } from '../../modules/local/medaka_variant'
include { TABIX_BGZIP as MEDAKA_BGZIP_VCF       } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_TABIX as MEDAKA_TABIX_VCF       } from '../../modules/nf-core/modules/tabix/tabix/main'
//include { DEEPVARIANT                           } from '../../modules/local/deepvariant'
//include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF  } from '../../modules/nf-core/modules/tabix/tabix/main'
//include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF } from '../../modules/nf-core/modules/tabix/tabix/main'
include { PEPPER_MARGIN_DEEPVARIANT             } from '../../modules/local/pepper_margin_deepvariant'

include { RUN_DEEPVARIANT                       } from '../../subworkflows/local/vc_deepvariant'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai
    ch_intervals
    ch_intervals_bed_gz_tbi
    ch_intervals_bed_combine_gz_tbi
    ch_intervals_bed_combine_gz

    main:
    ch_short_calls_vcf              = Channel.empty()
    ch_short_calls_vcf_tbi          = Channel.empty()
    ch_short_calls_gvcf             = Channel.empty()
    ch_short_calls_gvcf_tbi         = Channel.empty()
    ch_versions                     = Channel.empty()

    /*
     * Remap channels with intervals
     */
    view_sortbam_intervals = ch_view_sortbam.combine(ch_intervals)
        .map{ meta, input, index, intervals, num_intervals ->
            new_meta = meta.clone()
            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_meta.sample = meta.id
            new_meta.id = num_intervals <= 1 ? meta.id : meta.id + "_" + intervals.baseName
            new_meta.num_intervals = num_intervals

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [new_meta, input, index, intervals_new]
        }

    /*
     * Call short variants
     */
    if (params.variant_caller == 'medaka') {

        /*
         * Call short variants with medaka
         */
        MEDAKA_VARIANT( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(medaka_version = MEDAKA_VARIANT.out.versions)

        /*
         * Zip medaka vcf
         */
        MEDAKA_BGZIP_VCF( MEDAKA_VARIANT.out.vcf )
        ch_short_calls_vcf  = MEDAKA_BGZIP_VCF.out.gz
        ch_versions = ch_versions.mix(bgzip_version = MEDAKA_BGZIP_VCF.out.versions)

        /*
         * Index medaka vcf.gz
         */
        MEDAKA_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = MEDAKA_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(tabix_version = MEDAKA_TABIX_VCF.out.versions)

    } else if (params.variant_caller == 'deepvariant') {
        RUN_DEEPVARIANT(view_sortbam_intervals, ch_fasta, ch_fai, ch_intervals, ch_intervals_bed_gz_tbi, ch_intervals_bed_combine_gz_tbi, ch_intervals_bed_combine_gz)
        
        ch_short_calls_vcf = RUN_DEEPVARIANT.out.deepvariant_vcf
        ch_versions        = ch_versions.mix(RUN_DEEPVARIANT.out.versions)

    } else {

        /*
         * Call variants with pepper_margin_deepvariant (automatic zip + index, docker + singularity only)
         */
        PEPPER_MARGIN_DEEPVARIANT( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf = PEPPER_MARGIN_DEEPVARIANT.out.vcf
        ch_short_calls_vcf_tbi = PEPPER_MARGIN_DEEPVARIANT.out.tbi
        ch_versions = ch_versions.mix(PEPPER_MARGIN_DEEPVARIANT.out.versions)
    }

    emit:
    ch_short_calls_vcf
    ch_short_calls_vcf_tbi
    ch_short_calls_gvcf
    ch_short_calls_gvcf_tbi
    ch_versions
}
