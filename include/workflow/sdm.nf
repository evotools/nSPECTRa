include { get_individuals; get_breeds } from "../process/prerun"
include { makeAnnotation; annotateVcf } from '../process/annotate'
include { sdm; filter_sdm; count_sdm } from "../process/sdm"

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

        // Generate outputs
        count_sdm( filter_sdm.out )
}