include { get_individuals; get_breeds } from "../process/prerun"
include { makeAnnotation; annotateVcf } from '../process/annotate'
include { sdm; filter_sdm; count_sdm } from "../process/sdm"
include { make_ksfs; sdm_plot } from "../process/sdm"
include { repeat_mask_split_sdm} from "../process/sdm"

workflow SDM {
    take:
        vcf
        tbi
        ancfa
        ancfai
        reffasta
        reffai
        masks_ch
        chromosomeList

    main:

        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        // Get annotation
        makeAnnotation(vcf, ancfa, ancfai)
        annotateVcf(vcf, tbi, masks_ch, makeAnnotation.out[0], makeAnnotation.out[1])
        
        // Get individuals' ids
        // get_individuals( annotateVcf.out[0], annotateVcf.out[1] )

        // Import breeds' lists 
        breeds_ch = Channel
            .fromPath("${params.pops_folder}/*.txt")
            .map { file -> tuple(file.simpleName, file) }
        
        // Combine chromosomes and breeds
        combined_ch = breeds_ch.combine(chromosomes_ch)

        // Run dinuc pipeline
        sdm( annotateVcf.out[0], annotateVcf.out[1], reffasta, reffai, combined_ch )

        // Filter results
        filter_sdm(breeds_ch, sdm.out.collect() )

        // Get data in/out rm
        repeat_mask_split_sdm(filter_sdm.out.bed.combine(masks_ch))

        // Generate outputs
        count_sdm( filter_sdm.out.rdata )

        // Prepare Ksfs files
        make_ksfs(breeds_ch, sdm.out.collect())

        // Make plots for sdm results
        all_counts = count_sdm.out[0]
        //sdm_plot( breeds_ch, all_counts.collect() ) 
        sdm_plot( all_counts.collect() ) 
}