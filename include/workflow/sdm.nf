include { get_individuals; get_breeds; } from "../process/prerun"
include { sdm; filter_sdm; count_sdm } from "../process/sdm"
include { make_ksfs; sdm_plot } from "../process/sdm"
include { repeat_mask_split_sdm} from "../process/sdm"
include { sdm_matrix } from '../process/sdm.nf'

workflow SDM {
    take:
        vcf_by_chr
        reffasta
        reffai
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
        count_sdm( filtered_ch.rdata )

        // Create matrix of SDMs
        count_sdm.out[0] | collect | sdm_matrix

        // Prepare Ksfs files
        raw_sdm | make_ksfs

        // Make plots for sdm results
        all_counts = count_sdm.out[0]
        //sdm_plot( breeds_ch, all_counts.collect() ) 
        sdm_plot( all_counts.collect() ) 
}