
include { chromosomeList } from '../process/prerun'
include { mutyper; group_results } from '../process/mutyper'
include { mutyper_full; ksfs } from '../process/mutyper'
include { mutyper_variant; mutyper_spectra; mutyper_concat } from '../process/mutyper'
include { mutyper_full_parallel; plot_results } from '../process/mutyper'
include { count_mutations; count_mutations_csq } from '../process/mutyper'
include { kmercount; normalize_results } from '../process/mutyper'
include { extract_tags; extract_csq } from '../process/mutyper'
include { combine_counts; combine_csqs } from '../process/mutyper'


// Workflow
workflow MUTYPER {
    take:
        vcf
        tbi
        anc_fa
        anc_fai
        chromosomeList
        masks_ch
        vcf_by_chr
        // ne_time_ch

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
        combined_ch = vcf_by_chr.combine(k_list)

        // Run meryl counter
        kmercount( anc_fa, anc_fai, Channel.from(k_list) ) 

        // Run mutyper variant on each chromosome
        mutyper_variant(combined_ch, anc_fa.collect(), anc_fai.collect(), masks_ch.collect())
        // Then extract the spectrum
        | mutyper_spectra
        // In the meanwhile, group the mutyper VCFs and concatenate them
        mutyper_variant.out | groupTuple(by: 0) | mutyper_concat

        /* Generate mutyper-annotated vcf */
        // mutyper_full( vcf, tbi, anc_fa, anc_fai, masks_ch, Channel.from(k_list) )
        // mutyper_full_parallel( vcf, tbi, anc_fa, anc_fai, masks_ch, chromosomeList, Channel.from(k_list) )

        /* Generate the counts manually */
        // mutyper_full_parallel.out
        mutyper_variant.out
        | filter{ it[0].toInteger() < 8 }
        | extract_tags
        | map{
            k, region, tsv ->
            [k, region, tsv, file("${baseDir}/assets/K${k}_mutations.txt")]
        }
        | count_mutations
        | groupTuple(by: 0)
        | combine_counts

        // If VEP ran, compare changes by consequence
        if (params.vep){
            // mutyper_full_parallel.out
            mutyper_variant.out
            | filter{ it[0].toInteger() < 8 }
            | extract_csq
            | map{
                k, region, tsv ->
                [k, region, tsv, file("${baseDir}/assets/K${k}_mutations.txt"), file("${baseDir}/assets/VEPpriority")]
            }
            | count_mutations_csq
            | groupTuple(by: 0)
            | combine_csqs
        }

        /* Start mutyper on each chromosome separately */
        // mutyper( combined_ch, anc_fa.collect(), anc_fai.collect(), masks_ch.collect() )

        /* Collect outputs */
        // group_results( mutyper.out.groupTuple(by : 0) )
        group_results( mutyper_spectra.out.groupTuple(by : 0) )
        plot_results( group_results.out )

        // Normalize results
        // normalize_results( group_results.out, kmercount.out.collect() )

        // Grep list of samples, and drop smallest groups
        ksfs_inputs = Channel.fromPath("${params.pops_folder}/*.txt")
            | map { file -> tuple(file.simpleName, file, file.countLines() ) }
            | filter { it -> it[2] >= params.min_pop_size }
            | map { it -> tuple( it[0], it[1] ) }
            | combine( mutyper_concat.out )
            // | combine( mutyper_full_parallel.out )

        // Extract ksfs vector for each breed
        ksfs( ksfs_inputs )
        
    // emit:
    //     mutyper.out
    //     group_results.out
}