
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

    stub:
    """
    echo "seq1.r"
    echo "seq2.r"
    echo "seq3.r"
    """

    script:
    """
    tabix -l input.vcf.gz
    """
}

process filtering {
    tag "filter"
    label "medium"

    input:
    path variants
    path variants_idx 

    output:
    path "missingness_het_${params.reference}.txt"

    
    stub:
    """
    touch missingness_het_${params.reference}.txt
    """

    script:
    """
    if [[ ${variants} == *.bcf ]]; then
        plink --cow --bcf ${variants} --allow-extra-chr --mind 0.1 --double-id --het --out missingness_het_${params.reference} --threads ${task.cpus}
    else
        plink --cow --vcf ${variants} --allow-extra-chr --mind 0.1 --double-id --het --out missingness_het_${params.reference} --threads ${task.cpus}
    fi
    awk 'NR>1 && \$3/\$5>0.05 && \$3/\$5<0.95 {print \$1}' missingness_het_${params.reference}.het > missingness_het_${params.reference}.txt
    """
}

process keep_biallelic_snps {
    tag "bisnp"
    label "medium_smallmem_parallel"

    input:
    tuple path(variants), path(variants_idx), val(contig)

    output:
    tuple val(contig), path("biallelic_snps.${contig}.vcf.gz"), path("biallelic_snps.${contig}.vcf.gz.tbi")

    
    stub:
    """
    touch biallelic_snps.${contig}.vcf.gz
    touch biallelic_snps.${contig}.vcf.gz.tbi
    """

    script:
    """
    bcftools annotate -x INFO --threads ${task.cpus} ${variants} | \
        bcftools view --threads ${task.cpus} -t ${contig} -m 2 -M 2 -v snps -O z > biallelic_snps.${contig}.vcf.gz && tabix -p vcf biallelic_snps.${contig}.vcf.gz
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

    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """

    script:
    """
    bedtools intersect -header -u -a $variants -b ${extract} | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
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

    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """

    script:
    """
    bedtools intersect -header -v -a $variants -b ${exclude} | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
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

    
    stub:
    """
    touch prefilter.vcf.gz
    touch prefilter.vcf.gz.tbi
    """

    script:
    """
    bedtools intersect -v -a ${extract} -b ${exclude} > shortlisted.txt
    bedtools intersect -header -u -a $variants -b shortlisted.txt | bgzip -c > prefilter.vcf.gz && \
        tabix -p vcf prefilter.vcf.gz
    """
}


process apply_filter {
    tag "filter"
    label "medium_smallmem_parallel"

    input:
    tuple val(chrom), path(variants), path(variants_idx)
    path id_to_keep 

    output:
    tuple val(chrom), path("filtered.vcf.gz"), path("filtered.vcf.gz.tbi")
    

    
    stub:
    """
    touch filtered.${chrom}.vcf.gz
    touch filtered.${chrom}.vcf.gz.tbi
    """

    script:
    """
    bcftools view --threads ${task.cpus} -S ${id_to_keep} -m 2 -M 2 -q 0.05:minor -O z ${variants} > filtered.${chrom}.vcf.gz 
    tabix -p vcf filtered.${chrom}.vcf.gz
    """
}


process select_noncoding {
    label "medium_smallmem_parallel"

    input:
        tuple val(chrom), path(vcf), path(tbi)
        
    output:
        tuple val(chrom), path("${vcf.simpleName}.${chrom}.noncoding.vcf.gz"), path("${vcf.simpleName}.${chrom}.noncoding.vcf.gz.tbi")

    stub:
    """
    touch ${vcf.simpleName}.${chrom}.noncoding.vcf.gz
    touch ${vcf.simpleName}.${chrom}.noncoding.vcf.gz.tbi
    """

    script:
    """
    bcftools view --threads ${task.cpus} -O z -i 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.${chrom}.noncoding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.${chrom}.noncoding.vcf.gz
    """
}

process daf_filter {
    label "medium_smallmem_parallel"
    afterScript "rm tmp.vcf.gz*"

    input:
        tuple val(chrom), path('variants.vcf.gz'), path('variants.vcf.gz.tbi')

    output:
        tuple val(chrom), path("variants_DAF_${chrom}.vcf.gz"), path("variants_DAF_${chrom}.vcf.gz.tbi")

    stub:
    """
    touch variants_DAF_${chrom}.vcf.gz
    touch variants_DAF_${chrom}.vcf.gz.tbi
    """

    script:
    """
    bcftools +fill-tags variants.vcf.gz  -- -t AF,AC,AN |\
        bcftools filter --threads ${task.cpus} -O z -i 'INFO/AA != "-"' > tmp.vcf.gz && \
        bcftools index -t tmp.vcf.gz
    # Compute DAF
    ComputeDAF -v tmp.vcf.gz |\
        bcftools filter --threads ${task.cpus} -O z -i 'INFO/DAF <= ${params.max_derivate_allele_freq}' > variants_DAF_${chrom}.vcf.gz && \
        bcftools index --threads ${task.cpus} -t variants_DAF_${chrom}.vcf.gz
    """
}

process select_coding {
    label "medium_smallmem_parallel"

    input:
        tuple val(chrom), path(vcf), path(tbi)

    output:
        tuple val(chrom), path("${vcf.simpleName}.${chrom}.coding.vcf.gz"), path("${vcf.simpleName}.${chrom}.coding.vcf.gz.tbi")

    stub:
    """
    touch ${vcf.simpleName}.${chrom}.coding.vcf.gz
    touch ${vcf.simpleName}.${chrom}.coding.vcf.gz.tbi
    """

    script:
    """
    bcftools view --threads ${task.cpus} -O z -e 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.${chrom}.coding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.${chrom}.coding.vcf.gz
    """
}
