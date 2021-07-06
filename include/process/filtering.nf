
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