
/*
 * Phase 6: variants filtering
 */
process get_sequences {
    label "small"

    input:
    path "input.vcf.gz"
    path "input.vcf.gz.tbi"

    output:
    stdout

    script:
    """
    tabix -l input.vcf.gz
    """

    stub:
    """
    echo "seq1.r"
    echo "seq2.r"
    echo "seq3.r"
    """
}


process keep_biallelic_snps {
    tag "bisnp"
    label "medium_smallmem_parallel"

    input:
    tuple path(variants), path(variants_idx), val(contig)

    output:
    tuple val(contig), path("biallelic_snps.${contig}.vcf.gz"), path("biallelic_snps.${contig}.vcf.gz.tbi")

    script:
    """
    bcftools annotate -x INFO -r ${contig} --threads ${task.cpus} ${variants} | \
        bcftools view --threads ${task.cpus} -m 2 -M 2 -v snps -O z > biallelic_snps.${contig}.vcf.gz && tabix -p vcf biallelic_snps.${contig}.vcf.gz
    """
    
    stub:
    """
    touch biallelic_snps.${contig}.vcf.gz
    touch biallelic_snps.${contig}.vcf.gz.tbi
    """
}

process extract {
    tag "extract"
    label "medium"

    input:
    path variants
    path variants_idx 
    path extract

    output:
    path "prefilter.vcf.gz"
    path "prefilter.vcf.gz.tbi"

    script:
    """
    bedtools intersect -header -u -a $variants -b ${extract} | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
    """
    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """
}

process exclude {
    tag "exclude"
    label "medium"

    input:
    path variants
    path variants_idx 
    path exclude

    output:
    path "prefilter.vcf.gz"
    path "prefilter.vcf.gz.tbi"

    script:
    """
    bedtools intersect -header -v -a $variants -b ${exclude} | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
    """
    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """
}

process extract_exclude {
    tag "extexc"
    label "medium"

    input:
    path variants
    path variants_idx 
    path extract
    path exclude

    output:
    path "prefilter.vcf.gz"
    path "prefilter.vcf.gz.tbi"

    script:
    """
    bedtools intersect -v -a ${extract} -b ${exclude} > shortlisted.txt
    bedtools intersect -header -u -a $variants -b shortlisted.txt | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
    """
    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """
}


process select_noncoding {
    label "medium_smallmem_parallel"

    input:
        tuple val(chrom), path(vcf), path(tbi)
        
    output:
        tuple val(chrom), path("${vcf.simpleName}.${chrom}.noncoding.vcf.gz"), path("${vcf.simpleName}.${chrom}.noncoding.vcf.gz.tbi")

    script:
    """
    bcftools view --threads ${task.cpus} -O z -i 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.${chrom}.noncoding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.${chrom}.noncoding.vcf.gz
    """

    stub:
    """
    touch ${vcf.simpleName}.${chrom}.noncoding.vcf.gz
    touch ${vcf.simpleName}.${chrom}.noncoding.vcf.gz.tbi
    """
}

process daf_filter {
    label "medium_smallmem_parallel"

    input:
        tuple val(chrom), path('variants.vcf.gz'), path('variants.vcf.gz.tbi')

    output:
        tuple val(chrom), path("variants_DAF_${chrom}.vcf.gz"), path("variants_DAF_${chrom}.vcf.gz.tbi")

    script:
    // Filter alleles strictly, if required
    def filter_discordant = params.strict_allele_matching ? "-e 'REF!=AA && ALT!=AA'" : "" 
    if (params.max_derivate_allele_freq)
    """
    # Keep sites with AA, then with DAF<= max value
    bcftools filter --threads ${task.cpus} -O u -i 'AA != "-"' variants.vcf.gz | \
        bcftools filter --threads ${task.cpus} -O u ${filter_discordant} variants.vcf.gz | \
        bcftools filter --threads ${task.cpus} -O z -i 'INFO/DAF <= ${params.max_derivate_allele_freq}' > variants_DAF_${chrom}.vcf.gz && \
        bcftools index --threads ${task.cpus} -t variants_DAF_${chrom}.vcf.gz
    """
    else
    """
    # Keep sites with AA, then with DAF<= max value
    bcftools filter --threads ${task.cpus} -O u -i 'AA != "-"' variants.vcf.gz | \
        bcftools filter --threads ${task.cpus} -O z ${filter_discordant} > variants_DAF_${chrom}.vcf.gz && \
        bcftools index --threads ${task.cpus} -t variants_DAF_${chrom}.vcf.gz
    """

    stub:
    """
    touch variants_DAF_${chrom}.vcf.gz
    touch variants_DAF_${chrom}.vcf.gz.tbi
    """
}

process select_coding {
    label "medium_smallmem_parallel"

    input:
        tuple val(chrom), path(vcf), path(tbi)

    output:
        tuple val(chrom), path("${vcf.simpleName}.${chrom}.coding.vcf.gz"), path("${vcf.simpleName}.${chrom}.coding.vcf.gz.tbi")

    script:
    """
    bcftools view --threads ${task.cpus} -O z -e 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.${chrom}.coding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.${chrom}.coding.vcf.gz
    """

    stub:
    """
    touch ${vcf.simpleName}.${chrom}.coding.vcf.gz
    touch ${vcf.simpleName}.${chrom}.coding.vcf.gz.tbi
    """
}
