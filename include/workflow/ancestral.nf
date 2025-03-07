
include { hal2maf } from '../process/hal2maf'
include { rename_fasta as rename_ref; rename_fasta as rename_anc } from '../process/maf2fasta'
include {maf2bed; bed2vbed; bed2ancfa; collectAncfa; filter_by_size } from '../process/maf2fasta'
include { makeRefTgtFasta; splitfasta; makefai; makefai as makefai_ref } from '../process/maf2fasta'
include { liftToAncestor; makeAncestralSequence } from '../process/liftToAncestor'
include { makeAnnotation } from '../process/annotate'
include { chromosomeList } from '../process/prerun'
include { get_hal } from '../process/dependencies'

// Workflow
workflow ANCESTRAL {
    main:
        /*Create ancestral fasta file*/
        // hal to maf, needed for constrained elements

        if (params.ancestral_fna){
            // Import ancestral fasta
            anc_fa = Channel.fromPath(params.ancestral_fna)
            anc_fa | view{"Using ancestral genome: ${it}"}
            if (file("${params.ancestral_fna}.fai").exists()){
                anc_fai = Channel.fromPath("${params.ancestral_fna}.fai")
                anc_fai | view{"Using ancestral genome index: ${it}"}
            } else {
                anc_fai = makefai(anc_fa)
            }
            
            // Extract the different genomes and split it into chunks to speed up the process
            if (params.reference_fna){
                ch_ref = Channel.fromPath(params.reference_fna)
                ch_ref | view{"Using reference genome: ${it}"}
                if (file("${params.reference_fna}.fai").exists()){
                    ch_ref_fai = Channel.fromPath("${params.reference_fna}.fai")
                    ch_ref_fai | view{"Using reference genome index: ${it}"}
                } else {
                    ch_ref_fai = makefai_ref(ch_ref)
                }
                if (params.ref_min_size){
                    filter_by_size(ch_ref, ch_ref_fai)
                    ch_ref = filter_by_size.out[0]
                    ch_ref_fai = filter_by_size.out[1]
                }
            } else {
                if (params.hal) { ch_hal = file(params.hal) } else { exit 1, 'Hal file not specified!' }
                // hal2maf( ch_hal )
                makeRefTgtFasta( ch_hal )
                if (params.ref_min_size){
                    filter_by_size(makeRefTgtFasta.out[0], makeRefTgtFasta.out[2])
                    ch_ref = filter_by_size.out[0]
                    ch_ref_fai = filter_by_size.out[1]
                } else {
                    ch_ref = makeRefTgtFasta.out[0]
                    ch_ref_fai = makeRefTgtFasta.out[2]
                }
            }
        } else {
            if (params.hal) { ch_hal = file(params.hal) } else { exit 1, 'Hal file not specified!' }
            hal2maf( ch_hal )

            // maf to bed
            maf2bed(hal2maf.out)

            // Extract the different genomes and split it into chunks to speed up the process
            makeRefTgtFasta( ch_hal )
            if (params.ref_min_size){
                filter_by_size(makeRefTgtFasta.out[0], makeRefTgtFasta.out[2])
                ch_ref = filter_by_size.out[0]
                ch_ref_fai = filter_by_size.out[1]
            } else {
                ch_ref = makeRefTgtFasta.out[0]
                ch_ref_fai = makeRefTgtFasta.out[2]
            }

            sequences = ch_ref.splitFasta(record: [ id: true, header: true ], by: 1)
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
        anc_fna = anc_fa.collect()
        anc_fai = anc_fai.collect()
        ref_fna = ch_ref.collect()
        ref_fai = ch_ref_fai.collect()
}