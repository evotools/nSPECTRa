
include {
    hal4d;
    phyloFit;
    phyloPtrain;
    combine_bed;
    extract_conserved;
    vcf_drop_conserved;
    phyloP;
    make4dmaf;
    msa_view;
    hal_genomes;
    wig2bedgraph;
    halTree;
} from '../process/phylop'


workflow CONSTRAINED {
    take:
        vcf
        tbi
        chromosomeList

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

        // Get chromosome list
        chromosomes_ch = chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }

        // Extract 4d elements
        maf4d = hal4d(hal_ch, exons_ch)
        | combine(hal_genomes(hal_ch))
        | map{
            hal, maf, genomes_env ->
            String genomes = genomes_env as String
            [hal, maf, genomes]
        }
        | make4dmaf
        | msa_view

        // Fit model using 4D codons
        model_ch = maf4d
            | combine(hal_ch | halTree)
            | phyloFit

        // Run PhyloP
        phylop_ch = chromosomes_ch
        | combine(hal_ch)
        | combine(model_ch)
        | phyloP
        | wig2bedgraph
        // combine selected
        | collect
        | combine_bed

        // Extract conserved
        conserved_ch = phylop_ch | extract_conserved

        // Perform actual filtering
        vcf_drop_conserved(vcf, tbi, conserved_ch)

    emit:
        vcf_drop_conserved.out.vcf
        vcf_drop_conserved.out.tbi
}