
/*
 * Phase 6: variants filtering
 */
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
    label "small"

    input:
    path variants
    path variants_idx
    path id_to_keep 

    output:
    path "filtered.vcf.gz"
    path "filtered.vcf.gz.tbi"
    

    
    stub:
    """
    touch filtered.vcf.gz
    touch filtered.vcf.gz.tbi
    """

    script:
    """
    bcftools view -S ${id_to_keep} -m 2 -M 2 -q 0.05:minor -O z ${variants} > filtered.vcf.gz 
    tabix -p vcf filtered.vcf.gz
    """
}


process select_noncoding {
    tag "ncd"
    label "small"

    input:
        path vcf
        path tbi

    output:
        path "${vcf.simpleName}.noncoding.vcf.gz"
        path "${vcf.simpleName}.noncoding.vcf.gz.tbi"

    script:
    """
    bcftools view -O z -i 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.noncoding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.noncoding.vcf.gz
    """
}

process select_coding {
    tag "cd"
    label "small"

    input:
        path vcf
        path tbi

    output:
        path "${vcf.simpleName}.coding.vcf.gz"
        path "${vcf.simpleName}.coding.vcf.gz.tbi"

    script:
    """
    bcftools view -O z -e 'CSQ[*] ~ "intergenic_variant" || CSQ[*] ~ "intron_variant"' ${vcf} > ${vcf.simpleName}.coding.vcf.gz && \
        tabix -p vcf ${vcf.simpleName}.coding.vcf.gz
    """
}