
include { chromosomeList } from '../process/prerun'
include { mutyper; group_results; plot_results } from '../process/mutyper'
include { mutyper_full; ksfs } from '../process/mutyper'


// Workflow
workflow MUTYPER {
    take:
        vcf
        tbi
        anc_fa
        anc_fai
        chromosomeList
        masks_ch
        ne_time_ch

    main:
        // Get chromosome list
        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        // Create list of k
        k_list = params.k?.tokenize(',').flatten()

        // Define sequence ids
        // chromosomeList( vcf, tbi )
        combined_ch = chromosomes_ch.combine(k_list)

        /* Generate mutyper-annotated vcf */
        mutyper_full( vcf, tbi, anc_fa, anc_fai, masks_ch, Channel.from(k_list) )

        /* Start mutyper on each chromosome separately */
        mutyper( vcf, tbi, anc_fa, anc_fai, masks_ch, combined_ch )

        /* Collect outputs */
        group_results( mutyper.out.groupTuple(by : 0) )
        plot_results( group_results.out )

        // Grep list of samples, and drop smallest groups
        ksfs_inputs = Channel
            .fromPath("${params.pops_folder}/*.txt")
            .map { file -> tuple(file.simpleName, file, file.countLines() ) }
            .filter { it -> it[2] >= params.min_pop_size }
            .map { it -> tuple( it[0], it[1] ) }
            .combine( mutyper_full.out )

        // Extract ksfs vector for each breed
        ksfs( ksfs_inputs )
        
    // emit:
    //     mutyper.out
    //     group_results.out
}