
include { hal2maf } from '../process/hal2maf'
include {maf2bed; bed2vbed; bed2ancfa; collectAncfa } from '../process/maf2fasta'
include { makeRefTgtFasta; splitfasta; makefai } from '../process/maf2fasta'
include { liftToAncestor; makeAncestralSequence } from '../process/liftToAncestor'
include { makeAnnotation } from '../process/annotate'
include { chromosomeList } from '../process/prerun'
include { get_hal } from '../process/dependencies'

// Workflow
workflow ANCESTRAL {
    main:
        /*Create ancestral fasta file*/
        // Get hal
        get_hal()

        if (params.ancestral){
            // Import ancestral fasta
            anc_fa = file(params.ancestral)
            makefai(anc_fa)
            anc_fai = makefai.out
            
            // Extract the different genomes and split it into chunks to speed up the process
            makeRefTgtFasta( get_hal.out )
        } else {
            // hal to maf
            hal2maf( get_hal.out )
            // maf to bed
            maf2bed(hal2maf.out)

            // Extract the different genomes and split it into chunks to speed up the process
            makeRefTgtFasta( get_hal.out )
            splitfasta(makeRefTgtFasta.out[0], makeRefTgtFasta.out[2])
            //chunked_ref = maf2bed.out.combine(splitfasta.out)
            
            // Bed to vertical bed
            bed2vbed( maf2bed.out, splitfasta.out.flatten() )
            // Bed to ancestral fasta
            bed2ancfa( bed2vbed.out )
            // Collect ancestral fasta
            collectAncfa( bed2ancfa.out.collect() )
            anc_fa = collectAncfa.out[0]
            anc_fai = collectAncfa.out[1]
        }
        
    emit:
        anc_fa
        anc_fai
        makeRefTgtFasta.out[0]
        makeRefTgtFasta.out[2]
}