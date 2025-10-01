include { get_individuals; get_breeds; } from "../process/prerun"
include { sdm; filter_sdm; count_sdm } from "../process/sdm"
include { make_ksfs; sdm_plot } from "../process/sdm"
include { repeat_mask_split_sdm} from "../process/sdm"
include { fetch_sites } from '../process/sdm.nf'
include { sdm_matrix } from '../process/sdm.nf'
include { count_mutations } from '../process/mutyper.nf'

workflow SDM {
    take:
        vcf_by_chr
        reffasta
        reffai
        ancfasta
        ancfai
        masks_ch
        chromosomeList

    main:

        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        // Import breeds' lists 
        breeds_ch = Channel
            .fromPath("${params.pops_folder}/*.txt")
            .map { file -> tuple(file.simpleName, file) }
        
        // Prepare chunks
        combined_ch = breeds_ch
        | combine(
            vcf_by_chr
        )

        // Run dinuc pipeline
        raw_sdm = breeds_ch
        | combine(
            sdm( combined_ch, reffasta, reffai )
            | groupTuple(by: 0),
            by: 0
        )

        // Filter SDMs
        filtered_ch = raw_sdm | filter_sdm

        // Get data in/out rm
        repeat_mask_split_sdm(filtered_ch.bed.combine(masks_ch))

        // Generate outputs
        sdmcounts_ch = count_sdm( filtered_ch.rdata )

        // Create matrix of SDMs
        sdmcounts_ch.counts | collect | sdm_matrix

        // We collect first and second changes in a full list
        first_change = sdmcounts_ch.first_change
            | collectFile("${params.outdir}/sdm_first_change.txt")
        second_change = sdmcounts_ch.second_change
            | collectFile("${params.outdir}/sdm_second_change.txt")
        // Extract the changes
        fetch_sites( vcf_by_chr | collect(), ancfasta, ancfai, first_change | mix(second_change) )
        | map {
            vcf, tbi ->
            tuple( "3", vcf, tbi, file("${baseDir}/assets/K3_mutations.txt"), "sdm/${vcf.simpleName}" )
        }
        // Count individual mutations spectrums
        | count_mutations

        // Prepare Ksfs files
        raw_sdm | make_ksfs

        // Make plots for sdm results
        all_counts = sdmcounts_ch.counts
        //sdm_plot( breeds_ch, all_counts.collect() ) 
        sdm_plot( all_counts.collect() )
}