// workflow to generate the ibd windows
include {get_ref_ibd; get_merge_ibd} from '../process/dependencies'
include {ibd; make_map; merge_ibd} from '../process/ibd'

workflow IBD {
    take:
        vcf
        tbi
        chromosomeList

    main:
        // Get chromosome list
        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        // get refined-ibd if missing
        if (params.refinedibd){
            ch_refibd = file(params.refinedibd)
        } else {
            get_ref_ibd()
            ch_refibd = get_ref_ibd.out
        }

        // get merge-ibd if missing
        if (params.mergeibd){
            ch_mergeibd = file(params.mergeibd)
        } else {
            get_merge_ibd()
            ch_mergeibd = get_merge_ibd.out
        }

        // make map
        make_map(vcf, tbi)

        // run refined-ibd
        ibd(vcf, tbi, ch_refibd, chromosomes_ch)

        // Merge-ibd
        merge_ibd(ibd.out, vcf, tbi, ch_mergeibd, make_map.out)

    emit:
        merge_ibd.out

}