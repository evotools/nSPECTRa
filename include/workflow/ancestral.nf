
include { hal2maf } from '../process/hal2maf'
include { rename_fasta as rename_ref; rename_fasta as rename_anc } from '../process/maf2fasta'
include {maf2bed; bed2vbed; bed2ancfa; collectAncfa } from '../process/maf2fasta'
include { makeRefTgtFasta; splitfasta; makefai; makefai as makefai_ref } from '../process/maf2fasta'
include { liftToAncestor; makeAncestralSequence } from '../process/liftToAncestor'
include { makeAnnotation } from '../process/annotate'
include { chromosomeList } from '../process/prerun'
include { get_hal } from '../process/dependencies'

// Workflow
workflow ANCESTRAL {
    take:
        cactus
        
    main:
        /*Create ancestral fasta file*/
        // hal to maf, needed for constrained elements

        if (params.ancestral){
            // Import ancestral fasta
            anc_fa = file(params.ancestral)
            makefai(anc_fa)
            anc_fai = makefai.out
            
            // Extract the different genomes and split it into chunks to speed up the process
            if (params.ref_fasta){
                ch_ref = file(params.ref_fasta)
                makefai_ref(ch_ref)
                ch_ref_fai = makefai_ref.out
            } else {
                if (params.hal) { ch_hal = file(params.hal) } else { exit 1, 'Hal file not specified!' }
                hal2maf( ch_hal, cactus )
                makeRefTgtFasta( ch_hal, cactus )
                ch_ref = makeRefTgtFasta.out[0]
                ch_ref_fai = makeRefTgtFasta.out[2]
            }
        } else {
            if (params.hal) { ch_hal = file(params.hal) } else { exit 1, 'Hal file not specified!' }
            hal2maf( ch_hal, cactus )

            // maf to bed
            maf2bed(hal2maf.out)

            // Extract the different genomes and split it into chunks to speed up the process
            makeRefTgtFasta( ch_hal, cactus )
            ch_ref = makeRefTgtFasta.out[0]
            ch_ref_fai = makeRefTgtFasta.out[2]

            sequences = makeRefTgtFasta.out[0]
                                .splitFasta(record: [ id: true, header: true ], by: 1)
            // splitfasta(makeRefTgtFasta.out[0], makeRefTgtFasta.out[2])
            //chunked_ref = maf2bed.out.combine(splitfasta.out)
            
            // Bed to vertical bed
            // bed2vbed( maf2bed.out, splitfasta.out.flatten() )
            bed2vbed( maf2bed.out, makeRefTgtFasta.out[0], makeRefTgtFasta.out[2], sequences )
            // Bed to ancestral fasta
            bed2ancfa( bed2vbed.out )
            // Collect ancestral fasta
            collectAncfa( bed2ancfa.out.collect() )
            anc_fa = collectAncfa.out[0]
            anc_fai = collectAncfa.out[1]

            // Rename sequences after being extracted from the HAL 
            // archive to match names in the VCF/VEP cache
            if (params.rename_hal_sequences){
                conversion_table = file( params.rename_hal_sequences )
                rename_ref(ch_ref, conversion_table)
                rename_anc(anc_fa, conversion_table)
                ch_ref = rename_ref.out[0]
                ch_ref_fai = rename_ref.out[1]
                anc_fa = rename_anc.out[0]
                anc_fai = rename_anc.out[1]
            }

        }
        
    emit:
        anc_fa
        anc_fai
        ch_ref
        ch_ref_fai
}