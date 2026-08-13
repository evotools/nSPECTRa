
include {
    bgcFilter;
    catsort as catsort_tract;
    catsort as catsort_informative;
    combine_bed;
    create_maf;
    extract_conserved;
    GENOME_INTERVALS;
    hal4d;
    hal_genomes;
    halSize;
    halTree;
    make4dmaf;
    msa_view as maf_to_4d_ss;
    msa_view as maf_to_ss;
    phastBias;
    phyloFit;
    phyloP;
    phyloPtrain;
    vcf_drop_intervals;
    wig2bedgraph;
} from '../process/phylop'


workflow NEUTRAL_MODEL {
    main:
        // Load hal file
        if (params.hal) {
            if (file(params.hal).exists()){
                hal_ch = Channel.fromPath(params.hal)
            } else {
                exit 1, "Hal file ${params.hal} not found"
            }
        } else {
            exit 1, 'Hal file not specified!'
        }
        // load exon bed file
        if (file(params.exon_bed).exists()){
            exons_ch = Channel.fromPath(params.exon_bed)
        } else {
            exit 1, "Exon BED file ${params.exon_bed} not found"
        }
        // Hal genomes
        genomes_ch = hal_genomes(hal_ch)

        // Extract 4d elements
        maf4d = hal4d(hal_ch, exons_ch)
        | combine(genomes_ch)
        | map{
            hal, maf, genomes_env ->
            String genomes = genomes_env as String
            [hal, maf, genomes]
        }
        | make4dmaf
        | maf_to_4d_ss

        // Fit model using 4D codons
        model_ch = maf4d
            | combine(hal_ch | halTree)
            | phyloFit

        //
        // Make maf file for BGC calculation
        //
        // First, create the reference sizes file
        ref_sizes_ch = halSize(hal_ch)

        // then, we create the intervals and we parse them into a nested channel
        intervals_ch = GENOME_INTERVALS(ref_sizes_ch) | flatten
    emit:
        model = model_ch
        intervals = intervals_ch
        sizes = ref_sizes_ch
}


workflow CONSTRAINED {
    take:
        vcf_by_chr_ch
        model_ch
        intervals_ch
        ref_sizes_ch

    main:
        // Load hal file
        if (params.hal) {
            if (file(params.hal).exists()){
                hal_ch = Channel.fromPath(params.hal)
            } else {
                exit 1, "Hal file ${params.hal} not found"
            }
        } else {
            exit 1, 'Hal file not specified!'
        }
        // load exon bed file
        if (file(params.exon_bed).exists()){
            exons_ch = Channel.fromPath(params.exon_bed)
        } else {
            exit 1, "Exon BED file ${params.exon_bed} not found"
        }

        // Run PhyloP
        phylop_ch = phyloP(
            hal_ch | combine(intervals_ch),
            model_ch | collect
        )
        | combine(ref_sizes_ch)
        | wig2bedgraph
        // combine selected
        | collect

        // Extract conserved
        conserved_ch = combine_bed(phylop_ch, "phylop", "PHYLOP") | extract_conserved

        // Perform actual filtering
        vcf_out_ch = vcf_drop_intervals(vcf_by_chr_ch, conserved_ch, "non-conserved", "PHYLOP/")

    emit:
        vcf = vcf_out_ch
}


// workflow to compute the biased gene conversion (BGC) regions
workflow BGC {
    take:
        vcf_by_chr_ch
        model_ch
        intervals_ch
        ref_sizes_ch

    main:
        // Load hal file
        if (params.hal) {
            if (file(params.hal).exists()){
                hal_ch = Channel.fromPath(params.hal)
            } else {
                exit 1, "Hal file ${params.hal} not found"
            }
        } else {
            exit 1, 'Hal file not specified!'
        }
        // load exon bed file
        if (file(params.exon_bed).exists()){
            exons_ch = Channel.fromPath(params.exon_bed)
        } else {
            exit 1, "Exon BED file ${params.exon_bed} not found"
        }
        maf_ch = create_maf(hal_ch | combine(intervals_ch | flatten))

        // Run pastBias
        phastbias_ch = phastBias(maf_ch, model_ch | collect)

        // Collect bigwig outputs
        large_bed_ch = phastbias_ch.wig
        | combine(ref_sizes_ch)
        | wig2bedgraph
        | collect
        // Save tract regions
        tracts_bed_ch = catsort_tract (
            phastbias_ch.tracts_bed | collect,
            "bgc_regions.bed",
            "PHAST/gBGC/",
            "-k 1,1 -k2,2n"
        )
        // Save also informative regions
        informative_bed_ch = catsort_informative (
            phastbias_ch.informative_bed | collect,
            "bgc_informative.bed",
            "PHAST/gBGC/",
            "-k 1,1 -k2,2n"
        )

        // Combine all BGC regions
        bgc_ch = combine_bed(large_bed_ch, "bgc", "PHASTBIAS")

        // Perform actual filtering
        vcf_out_ch = vcf_drop_intervals(vcf_by_chr_ch, tracts_bed_ch, 'no-bgc', "PHASTBIAS/")

    emit:
        vcf = vcf_out_ch
}