
include {hal4d; phyloFit; phyloPtrain; combine_bed } from '../process/phylop'
include {phyloP; bigWigToBed; filter; phastCons} from '../process/phylop'
include {make4dmaf; msa_view; makeMaf; collect; combine_mafs} from '../process/phylop'
include {get_phast; get_hal} from '../process/dependencies' 

workflow CONSTRAINED {
    take:
        hal
        vcf
        tbi
        chromosomeList

    main:
        // Get chromosome list
        chromosomeList
            .splitCsv(header: ['N','chrom'])
            .map{ row-> tuple(row.N, row.chrom) }
            .set{ chromosomes_ch }

        makeMaf(hal, chromosomes_ch)
        if (params.hal4d && params.exon_bed){
            hal4d(hal)
            make4dmaf( hal4d.out, hal )
            msa_view( make4dmaf.out, hal )
            //phyloPtrain(hal, hal4d.out)
            //model_ch = phyloPtrain.out
            phyloFit(msa_view.out, hal)
            model_ch = phyloFit.out
        } else {
            //makeMaf(hal)
            combine_mafs(makeMaf.out.collect())
            phyloFit(combine_mafs.out, hal)
            model_ch = phyloFit.out
        }
        // Run PhyloP
        //phyloP(hal, model_ch, chromosomes_ch)
        phastCons(hal, makeMaf.out, model_ch)

        // Collect bed
        bedfiles = phastCons.out[0]
        collect( bedfiles.collect() )

        // Convert single bigwig to bed
        //bigWigToBed(phyloP.out)

        // Combine single-chromosome beds
        //combine_bed( phyloP.out.collect() )

        // Perform actual filtering
        filter(vcf, tbi, collect.out)

    emit:
        //vcf
        //tbi
        filter.out[0]
        filter.out[1]

}