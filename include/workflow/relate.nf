include { chromosomeList } from "../process/prerun"
include { relate_format; relate; relate_mut; relate_mut_chr } from '../process/relate'
include { relate_mut_finalise; relate_mut_chr_finalise; relate_avg_mut } from '../process/relate'
include { relate_mut_chr_pop; relate_chr_pop_mut_finalise } from '../process/relate'
include { make_mask; makepopfile } from '../process/relate'
include { relate_ne; make_relate_map; relate_plot_pop } from '../process/relate'
include { makeAnnotation } from '../process/annotate'
include { makefai } from '../process/maf2fasta'

// Relate main workflow
workflow RELATE {
    take:
        vcf
        tbi
        anc_fa
        anc_fai
        ref_fa
        ref_fai
        chromosomeList
        bedmask

    main:
        ch_relate = file(params.relate_path)

        // check if there is a popfile, otherwise create one
        if (params.poplabels) {
            popfile_ch = Channel.fromPath(params.poplabels)
        } else {
            makepopfile( file(params.pops_folder) )
            popfile_ch = makepopfile.out
        }


        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        // Grep list of samples, and drop smallest groups
        breeds_ch = Channel
            .fromPath("${params.pops_folder}/*.txt")
            .map { file -> tuple(file.simpleName, file ) }
        
        // Combine chromosomes and breeds
        combined_ch = breeds_ch.combine(chromosomes_ch)
        
        // Prepare mask fasta
        make_mask( anc_fa, bedmask )
        makefai( make_mask.out )

        // Prepare data for relate
        // relate_format(vcf, tbi, breeds_ch, chromosomes_ch)
        relate_format(vcf, tbi, anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch, ch_relate, chromosomes_ch)

        // Create relate map files
        make_relate_map( relate_format.out )

        // Run relate
        relate( make_relate_map.out, ch_relate )

        // Run relate Ne
        relate_ne( relate.out[0].collect(), relate.out[1].collect(), chromosomeList, popfile_ch, ch_relate )
        ne_out = relate_ne.out

        // Estimate mutation spectra
        //relate_mut( chromosomeList, relate_ne.out[0].collect(), anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch )
        
        // Estimate mutation spectra (single chr)
        relate_mut_chr( chromosomes_ch, ne_out[0].collect(), anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch, ch_relate )
        //relate_mut_chr( chromosomes_ch, relate_ne.out[0].collect(), anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch )

        // Finalize per-contig mutation rates
        relate_mut_chr_finalise( relate_mut_chr.out, chromosomeList, popfile_ch, ch_relate )

        // Collect finalised values
        relate_mut_finalise( relate_mut_chr_finalise.out.collect(), chromosomeList, popfile_ch, ch_relate )

        // Run relate mutation by breed by chr
        //relate_mut_chr_pop( combined_ch, relate_ne.out[0].collect(), anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch )
        relate_mut_chr_pop( combined_ch, ne_out[0].collect(), anc_fa, anc_fai, make_mask.out, makefai.out, popfile_ch, ch_relate )

        // Finalize per-contig mutation rates
        relate_chr_pop_mut_finalise( relate_mut_chr_pop.out.groupTuple(by: 0), chromosomeList, ch_relate )

        // Calculate average mutation rate for population
        //relate_avg_mut( relate_ne.out[0].collect(), chromosomeList, popfile_ch )
        relate_avg_mut( ne_out[0].collect(), chromosomeList, popfile_ch, ch_relate )

        // Plot the relate results
        // relate_plot_pop( relate_chr_pop_mut_finalise.out.collect(), popfile_ch )
        
    emit:
        //relate_ne.out[2] 
        ne_out[2] 
        relate_avg_mut.out


}