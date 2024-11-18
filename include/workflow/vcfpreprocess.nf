
include { extract; exclude } from '../process/filtering'
include { extract_exclude; select_noncoding; select_coding } from '../process/filtering'
include { keep_biallelic_snps; daf_filter; get_sequences } from '../process/filtering'
// include { filtering; apply_filter } from '../process/filtering'
include { vep; makeAnnotation; annotateVcf } from '../process/annotate'
include { shapeit5; shapeit4; beagle} from '../process/prerun' 
include { chromosomeList; combineVcf } from '../process/prerun'
include { daf; smile; chunking } from '../process/prerun'
include { get_beagle; get_vep_cache } from '../process/dependencies'


// Run pre-process workflow
workflow PREPROCESS {
    take:
        ch_var
        ch_var_idx
        ch_ref
        ch_ref_fai
        ch_anc
        ch_anc_fai
        masks_ch
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
            // } else if (params.imputation == 'shapeit5'){
                // phased_vcf = shapeit5( biallelic_ch )
            } else {
                phased_vcf = shapeit4( biallelic_ch )
            }
        }

        // Annotate the output VCFs if either vep or sdm are requested
        if (params.vep || params.sdm){
            annotated = vep( phased_vcf, ch_ref, ch_ref_fai, annot_ch )
        } else {
            annotated = phased_vcf
        }

        // Annotate the output VCFs
        if (params.coding && !params.noncoding){
            processed_ch = select_coding( annotated )
        } else if (params.noncoding && !params.coding) {
            processed_ch = select_noncoding( annotated )
        } else {
            processed_ch = annotated
        }

        // Add ancestral allele information
        makeAnnotation(processed_ch, ch_anc.collect(), ch_anc_fai.collect())
        processed_ch = annotateVcf(makeAnnotation.out, masks_ch.collect())

        // Generate derived allele frequency plots prior DAF filtering.
        processed_ch | daf | collect | smile

        // Filter by derived allele freq.
        if (params.max_derivate_allele_freq){
            processed_ch = processed_ch | daf_filter
        }

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

        // Prepare chunks
        vcfs_ch = processed_ch | multiMap{
            chrom, vcf, tbi ->
            vcfs: vcf
            tbis: tbi
        }
        chunks_ch = chunking(vcfs_ch.vcfs, vcfs_ch.tbis)
        | splitCsv(header: false, sep: '\t')
        | combine(processed_ch, by: 0)

    emit:
        vcf = vcf_ch
        tbi = tbi_ch
        chroms = chr_list
        vcf_by_chr = processed_ch
        chunks_ch = chunks_ch
}