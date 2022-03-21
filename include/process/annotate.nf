

/*
 * Phase 5: make annotation for variants 
 */
process makeAnnotation {
    tag "makeAnnot"
    label 'medium'

    input:
    path vcf
    path ancestral
    path ancestral_tbi

    output:
    path "Ancestral_annotation_${params.reference}.txt.gz", emit: annotation
    path "Ancestral_annotation_${params.reference}.txt.gz.tbi", emit: annottbi


    stub:
    """
    touch Ancestral_annotation_${params.reference}.txt.gz
    touch Ancestral_annotation_${params.reference}.txt.gz.tbi
    """

    script:
    """
    bedtools getfasta -fi ${ancestral} -bed ${vcf} -bedOut | \
        awk 'BEGIN{OFS="\t"};{print \$1,\$2,\$NF}' | \
        bgzip -c > Ancestral_annotation_${params.reference}.txt.gz && \
        tabix -b 2 -e 2 -s 1 Ancestral_annotation_${params.reference}.txt.gz
    """
}

// Annotate VCF
process annotateVcf {
    label "medium"

    input:
    path vcf
    path tbi
    path masks_ch
    path annotation
    path annotation_idx

    output:
    path "genotypes.with_ancestral_allele.vcf.gz"
    path "genotypes.with_ancestral_allele.vcf.gz.tbi"


    stub:
    """
    touch genotypes.with_ancestral_allele.vcf.gz
    touch genotypes.with_ancestral_allele.vcf.gz.tbi
    """

    script:
    """
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > header_${params.reference}.txt
    bedtools intersect -header -v -a ${vcf} -b ${masks_ch} | \
        bcftools annotate -a ${annotation} -h header_${params.reference}.txt -O z -c CHROM,POS,AA > genotypes.with_ancestral_allele.vcf.gz && \
        tabix -p vcf genotypes.with_ancestral_allele.vcf.gz

    """
}


// Make VEP process
process vep {
    tag "vep"
    label "vep"

    conda (params.enable_conda ? "bioconda::ensembl-vep=${params.vep_release_maj}" : null)

    input:
    tuple val(chrom), file(vcf)
    path reffasta
    path reffai
    path input_annot

    output:
    path "genotypes.${chrom}.vep.vcf.gz"
    path "genotypes.${chrom}.vep.vcf.gz.tbi"
    
    
    stub:
    """
    touch genotypes.${chrom}.vep.vcf.gz
    touch genotypes.${chrom}.vep.vcf.gz.tbi
    """

    script:
    if (params.custom_vep)
    """
    if [ ! -e ${input_annot.simpleName}.tbi ]; then tabix -p gff ${input_annot}; fi
    vep -i ${vcf} -o stdout --vcf --fork ${task.cpus} --gff ${input_annot} --fasta ${reffasta} --variant_class --sift b --nearest symbol --distance 200 | bgzip -c > genotypes.${chrom}.vep.vcf.gz
    tabix -p vcf genotypes.${chrom}.vep.vcf.gz
    """
    else
    """
    vep -i ${vcf} -o stdout --vcf --fork ${task.cpus} --species ${params.species} --variant_class --sift b --nearest symbol --distance 200 --offline --dir_cache ${input_annot} | bgzip -c > genotypes.${chrom}.vep.vcf.gz
    tabix -p vcf genotypes.${chrom}.vep.vcf.gz
    """
}