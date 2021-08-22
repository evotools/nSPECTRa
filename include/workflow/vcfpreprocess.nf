
include { filtering; apply_filter; extract; exclude; extract_exclude } from '../process/filtering'
include { vep } from '../process/annotate'
include { shapeit4; beagle; chromosomeList; combine } from '../process/prerun'
include { get_beagle; get_vep_cache } from '../process/dependencies'


// Run pre-process workflow
workflow PREPROCESS {
    take:
        ch_var
        ch_var_idx
        ch_ref
        ch_ref_fai
    main:
        // Split variants in chromosomes and impute them
        if (params.chr_list){
            Channel
                .fromPath(params.chr_list)
                .splitCsv(header: ['N','chrom'])
                .map{ row-> tuple(row.N, row.chrom) }
                .set{ chromosomes_ch }
            chr_list = file(params.chr_list)
        } else {
            chromosomeList( ch_var, ch_var_idx )
            chromosomeList.out
                .splitCsv(header: ['N','chrom'])
                .map{ row-> tuple(row.N, row.chrom) }
                .set{ chromosomes_ch }
            chr_list = chromosomeList.out
        }

        // Check vep cache
        if (params.download_cache) {
            annot_ch = get_vep_cache()
        } else if (params.vep_cache) { 
            annot_ch = file(params.vep_cache)
        } else if (params.gff) { 
            annot_ch = file(params.gff) 
        } else { 
            annot_ch = get_vep_cache() 
        }

        // Impute with beagle or shapeit4
        if (params.imputation == 'beagle'){
            if (params.beagle){ 
                beagle_ch = file(params.beagle) 
            } else {
                get_beagle()
                beagle_ch = get_beagle.out
            }
            beagle( chromosomes_ch, ch_var, ch_var_idx, beagle_ch )
            vep( beagle.out, ch_ref, ch_ref_fai, annot_ch )
        } else {
            shapeit4( chromosomes_ch, ch_var, ch_var_idx )
            vep( shapeit4.out, ch_ref, ch_ref_fai, annot_ch )
        }
        // Annotate the output VCFs
        // vep( vcf_ch, ch_ref, ch_ref_fai, vep_cache )
        vcf_ch = vep.out[0]
        tbi_ch = vep.out[1]

        // Combine results
        combine( vcf_ch.collect(), tbi_ch.collect() )
        vcf_ch = combine.out[0]
        tbi_ch = combine.out[1]

        // Filter if requested
        if (params.filter == 'y'){
            filtering(vcf_ch, tbi_ch)
            apply_filter(vcf_ch, tbi_ch, filtering.out)
            vcf_ch = apply_filter.out[0]
            tbi_ch = apply_filter.out[1]
        } 

        // Narrow to given regions
        if (params.exclude && !params.extract){
            ch_exc = file(params.exclude)
            exclude(vcf_ch, tbi_ch, ch_exc)
            vcf_ch = extract_exclude.out[0]
            tbi_ch = extract_exclude.out[1]
        } else if (!params.exclude && params.extract){
            ch_ext = file(params.extract)
            extract(vcf_ch, tbi_ch, ch_ext)
            vcf_ch = extract_exclude.out[0]
            tbi_ch = extract_exclude.out[1]
        } else if (params.exclude && params.extract){
            ch_exc = file(params.exclude)
            ch_ext = file(params.extract)
            extract_exclude(vcf_ch, tbi_ch, ch_ext, ch_exc)
            vcf_ch = extract_exclude.out[0]
            tbi_ch = extract_exclude.out[1]
        }

    emit:
        vcf_ch
        tbi_ch
        chr_list
}