
include { chromosomeList } from '../process/prerun'
include { group_results; plot_results; ksfs } from '../process/mutyper'
include { mutyper_variant; mutyper_spectra; mutyper_concat } from '../process/mutyper'
include { count_mutations; count_mutations_csq; consequence_table } from '../process/mutyper'
include { kmercount; normalize_results } from '../process/mutyper'


// Workflow
workflow MUTYPER {
    take:
        anc_fa
        anc_fai
        masks_ch
        vcf_by_chr

    main:
        // check if there is a popfile, otherwise create one
        breeds_file = Channel.fromPath("${params.pops_folder}/*.txt")
            | map {
                fname ->
                [fname.baseName, fname]
            }
            | splitCsv
            | map {
                pop, iid -> [pop] + iid
            }
            | collectFile {
                pop, sample ->
                [ "meta.txt", "${pop}\t${sample}\n" ]
            }

        // Create list of k
        String ks = params.k as String
        if (ks.contains(',')){
                k_list = Channel.of(ks.tokenize(',')).flatten()
        } else {
                k_list = Channel.of(ks)
        }

        // Define sequence ids
        combined_ch = vcf_by_chr.combine(k_list)

        // Run meryl counter
        kmercount( anc_fa, anc_fai, k_list ) 

        // Run mutyper variant on each chromosome
        mutyper_variant(combined_ch, anc_fa.collect(), anc_fai.collect(), masks_ch.collect())
        | mutyper_spectra  // Then extract the spectrum

        // In the meanwhile, group the mutyper VCFs and concatenate them
        mutyper_vcf_ch = mutyper_variant.out | groupTuple(by: 0) | mutyper_concat

        /* Generate the counts manually */
        mutyper_vcf_ch
        | filter{ it[0].toInteger() < 8 }
        | map{
            k, vcf_fn, tbi_fn ->
            [k, vcf_fn, tbi_fn, file("${baseDir}/assets/K${k}_mutations.txt")]
        }
        | count_mutations

        // If VEP ran, compare changes by consequence
        if (params.vep){
            // Create consequence table
            mutyper_vcf_ch | consequence_table
            // Create change-by-csq counts 
            mutyper_vcf_ch
            | filter{ it[0].toInteger() < 8 }
            | map{
                k, vcf_fn, tbi_fn ->
                [k, vcf_fn, tbi_fn, file("${baseDir}/assets/K${k}_mutations.txt"), file("${baseDir}/assets/VEPpriority")]
            }
            | count_mutations_csq
        }

        /* Collect outputs */
        group_results( mutyper_spectra.out.groupTuple(by : 0) )
        plot_results( group_results.out, breeds_file.collect() )

        // Normalize results
        normalize_results( count_mutations.out.combine(kmercount.out, by: 0) )

        // Grep list of samples, and drop smallest groups
        ksfs_inputs = Channel.fromPath("${params.pops_folder}/*.txt")
            | map { file -> tuple(file.simpleName, file, file.countLines() ) }
            | filter { it -> it[2] >= params.min_pop_size }
            | map { it -> tuple( it[0], it[1] ) }
            | combine( mutyper_concat.out )

        // Extract ksfs vector for each breed
        ksfs( ksfs_inputs )
        
}