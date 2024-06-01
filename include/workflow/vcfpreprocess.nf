
include { extract; exclude } from '../process/filtering'
include { extract_exclude; select_noncoding; select_coding } from '../process/filtering'
include { keep_biallelic_snps; daf_filter; get_sequences } from '../process/filtering'
// include { filtering; apply_filter } from '../process/filtering'
include { vep } from '../process/annotate'
include { shapeit4; beagle} from '../process/prerun' 
include { chromosomeList; combineVcf } from '../process/prerun'
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
        // Get sequences to process
        sequences = get_sequences(ch_var, ch_var_idx) | splitText | map { it.trim() }

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

        // Faster processing of the raw VCF, running by contig
        biallelic_ch = ch_var | combine(ch_var_idx) | combine(sequences) | keep_biallelic_snps
        if (params.imputed){
            phased_vcf = biallelic_ch
        } else {
            // Impute with beagle or shapeit4
            if (params.imputation == 'beagle'){
                if (params.beagle){ 
                    beagle_ch = file(params.beagle) 
                } else {
                    get_beagle()
                    beagle_ch = get_beagle.out
                }
                phased_vcf = beagle( biallelic_ch, beagle_ch )
            } else {
                phased_vcf = shapeit4( biallelic_ch )
            }
        }

        // Annotate the output VCFs
        vep( phased_vcf, ch_ref, ch_ref_fai, annot_ch )

        // Annotate the output VCFs
        if (params.coding && !params.noncoding){
            processed_ch = select_coding( vep.out )
        } else if (params.noncoding && !params.coding) {
            processed_ch = select_noncoding( vep.out )
        } else {
            processed_ch = vep.out
        }

        // // Filter if requested
        // if (params.filter){
        //     filtering(ch_var, ch_var_idx)
        //     processed_ch = apply_filter(processed_ch.combine(filtering.out))
        // }

        // Combine results
        processed_ch | map{ chrom, vcf, tbi -> [vcf, tbi]} | flatten | collect | combineVcf
        vcf_ch = combineVcf.out[0]
        tbi_ch = combineVcf.out[1]

        // Narrow to given regions
        if (params.exclude && !params.extract){
            ch_exc = Channel.fromPath(params.exclude)
            exclude(vcf_ch, tbi_ch, ch_exc)
            vcf_ch = exclude.out[0]
            tbi_ch = exclude.out[1]
        } else if (!params.exclude && params.extract){
            ch_ext = Channel.fromPath(params.extract)
            extract(vcf_ch, tbi_ch, ch_ext)
            vcf_ch = extract.out[0]
            tbi_ch = extract.out[1]
        } else if (params.exclude && params.extract){
            ch_exc = Channel.fromPath(params.exclude)
            ch_ext = Channel.fromPath(params.extract)
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