
include { filtering; apply_filter; extract; exclude} from '../process/filtering'
include { extract_exclude; select_noncoding; select_coding } from '../process/filtering'
include { keep_biallelic_snps } from '../process/filtering'
include { vep } from '../process/annotate'
include { shapeit4; beagle; split_vcf} from '../process/prerun' 
include { chromosomeList; combine } from '../process/prerun'
include { daf; smile } from '../process/prerun'
include { get_beagle; get_vep_cache } from '../process/dependencies'


// Run pre-process workflow
workflow PREPROCESS {
    take:
        ch_var
        ch_var_idx
        ch_ref
        ch_ref_fai
    main:
        // Keep biallelic SNPs only
        keep_biallelic_snps(ch_var, ch_var_idx)
        ch_var = keep_biallelic_snps.out[0]
        ch_var_idx = keep_biallelic_snps.out[1]

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

        if (params.imputed){
            split_vcf( chromosomes_ch, ch_var, ch_var_idx )
            phased_vcf = split_vcf.out
        } else {
            // Impute with beagle or shapeit4
            if (params.imputation == 'beagle'){
                if (params.beagle){ 
                    beagle_ch = file(params.beagle) 
                } else {
                    get_beagle()
                    beagle_ch = get_beagle.out
                }
                beagle( chromosomes_ch, ch_var, ch_var_idx, beagle_ch )
                phased_vcf = beagle.out
            } else {
                shapeit4( chromosomes_ch, ch_var, ch_var_idx )
                phased_vcf = shapeit4.out
            }
        }

        // Annotate the output VCFs
        vep( phased_vcf, ch_ref, ch_ref_fai, annot_ch )
        vcf_ch = vep.out[0]
        tbi_ch = vep.out[1]

        // Combine results
        combine( vcf_ch.collect(), tbi_ch.collect() )

        // Annotate the output VCFs
        if (params.coding && !params.noncoding){
            select_coding( combine.out[0], combine.out[1] )
            vcf_ch = select_coding.out[0]
            vcf_ch = select_coding.out[1]
        } else if (params.noncoding && !params.coding) {
            select_noncoding( combine.out[0], combine.out[1] )
            vcf_ch = select_noncoding.out[0]
            vcf_ch = select_noncoding.out[1]
        } else {
            vcf_ch = combine.out[0]
            tbi_ch = combine.out[1]
        }

        // Filter if requested
        if (params.filter){
            filtering(vcf_ch, tbi_ch)
            apply_filter(vcf_ch, tbi_ch, filtering.out)
            vcf_ch = apply_filter.out[0]
            tbi_ch = apply_filter.out[1]
        } 

        // Narrow to given regions
        if (params.exclude && !params.extract){
            ch_exc = file(params.exclude)
            exclude(vcf_ch, tbi_ch, ch_exc)
            vcf_ch = exclude.out[0]
            tbi_ch = exclude.out[1]
        } else if (!params.exclude && params.extract){
            ch_ext = file(params.extract)
            extract(vcf_ch, tbi_ch, ch_ext)
            vcf_ch = extract.out[0]
            tbi_ch = extract.out[1]
        } else if (params.exclude && params.extract){
            ch_exc = file(params.exclude)
            ch_ext = file(params.extract)
            extract_exclude(vcf_ch, tbi_ch, ch_ext, ch_exc)
            vcf_ch = extract_exclude.out[0]
            tbi_ch = extract_exclude.out[1]
        }

        daf(vcf_ch, tbi_ch) | smile

        // Filter by derived allele freq.
        if (params.max_derivate_allele_freq){
            daf_filter(vcf_ch, tbi_ch)
            vcf_ch = daf_filter.out[0]
            tbi_ch = daf_filter.out[1]
        }

    emit:
        vcf_ch
        tbi_ch
        chr_list
}