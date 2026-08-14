// workflow to generate the ibd windows
include {get_ref_ibd; get_merge_ibd} from '../process/dependencies'
include {ibd; make_map; merge_ibd} from '../process/ibd'

workflow IBD {
    take:
        vcf_by_chr

    main:
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

        // run refined-ibd
        merge_ibd(
                ibd(
                make_map(vcf_by_chr), 
                ch_refibd
            ),
            ch_mergeibd
        )

    emit:
        merge_ibd.out

}