
include {gone_run; gone_inputs; collectNe} from "../process/gone"
include {gone_get} from '../process/dependencies'

workflow GONE {
    take:
        vcf_ch
        tbi_ch
        
    main:
        // get GONE
        gone_get()

        // import list of pops
        populations_ch = Channel
            .fromPath("${params.pops_folder}/*.txt")
            .map { file -> tuple(file.simpleName, file, file.countLines() ) }
            .filter { it -> it[2] >= params.min_pop_size }
            .map { it -> tuple( it[0], it[1] ) }

        // prepare gone
        gone_inputs( vcf_ch, tbi_ch, populations_ch )

        // run Ne
        gone_run( gone_inputs.out, gone_get.out )

        // Collect Ne values
        collectNe( gone_run.out.collect() )

    emit:
        collectNe.out
}