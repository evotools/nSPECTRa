

/*
 * Phase 5: make annotation for variants 
 */
process makeAnnotation {
    tag "makeAnnot"
    label 'medium'

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path ancestral
    path ancestral_tbi

    output:
    tuple val(chrom),
        path(vcf),
        path(tbi),
        path("Ancestral_annotation_${chrom}.txt.gz"),
        path("Ancestral_annotation_${chrom}.txt.gz.tbi")

    stub:
    """
    touch Ancestral_annotation_${params.reference}.txt.gz
    touch Ancestral_annotation_${params.reference}.txt.gz.tbi
    """

    script:
    """
    bedtools getfasta -fi ${ancestral} -bed ${vcf} -bedOut | \
        awk 'BEGIN{OFS="\t"};{print \$1,\$2,\$NF}' | \
        bgzip -c > Ancestral_annotation_${chrom}.txt.gz && \
        tabix -b 2 -e 2 -s 1 Ancestral_annotation_${chrom}.txt.gz
    """
}

// Annotate VCF
process annotateVcf {
    label "medium"

    input:
    tuple val(chrom),
        path('input.vcf.gz'),
        path('input.vcf.gz.tbi'),
        path("Ancestral_annotation.txt.gz"),
        path("Ancestral_annotation.txt.gz.tbi")
    path masks_ch

    output:
    tuple val(chrom),
        path("genotypes.with_ancestral_allele_${chrom}.vcf.gz"),
        path("genotypes.with_ancestral_allele_${chrom}.vcf.gz.tbi")


    stub:
    """
    touch genotypes.with_ancestral_allele_${chrom}.vcf.gz
    touch genotypes.with_ancestral_allele_${chrom}.vcf.gz.tbi
    """

    script:
    """
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > header_${params.reference}.txt
    bcftools annotate -a Ancestral_annotation.txt.gz -O u -h header_${params.reference}.txt -O z -c CHROM,POS,AA input.vcf.gz | \
        bcftools view -T ^${masks_ch} -O z > genotypes.with_ancestral_allele_${chrom}.vcf.gz && \
        tabix -p vcf genotypes.with_ancestral_allele_${chrom}.vcf.gz
    """
}


// Make VEP process
process vep {
    container params.vep_cache_version != 'latest' && params.vep_cache_version ? "ensemblorg/ensembl-vep:release_${params}" : "ensemblorg/ensembl-vep:latest"
    conda {params.enable_conda ? "bioconda::ensembl-vep=${params.vep_release_maj}" : null}

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path reffasta
    path reffai
    path input_annot

    output:
    tuple val(chrom), path("genotypes.${chrom}.vep.vcf.gz"), path("genotypes.${chrom}.vep.vcf.gz.tbi")
    
    stub:
    """
    touch genotypes.${chrom}.vep.vcf.gz
    touch genotypes.${chrom}.vep.vcf.gz.tbi
    """

    script:
    def cachev = params.vep_cache_version ? "--cache_version ${params.vep_cache_version}" : ""
    def sift = params.no_sift ? "" : "--sift b"
    if (params.custom_vep)
    """
    if [ ! -e ${input_annot.simpleName}.tbi ]; then tabix -p gff ${input_annot}; fi
    vep -i ${vcf} -o stdout --vcf --fork ${task.cpus} --gff ${input_annot} --fasta ${reffasta} --variant_class ${sift} --nearest symbol --distance 200 | bgzip -c > genotypes.${chrom}.vep.vcf.gz
    tabix -p vcf genotypes.${chrom}.vep.vcf.gz
    """
    else
    """
    vep -i ${vcf} -o stdout --vcf --fork ${task.cpus} --species ${params.species} --variant_class ${sift} --nearest symbol --distance 200 --offline --dir_cache ${input_annot} ${cachev} | bgzip -c > genotypes.${chrom}.vep.vcf.gz
    tabix -p vcf genotypes.${chrom}.vep.vcf.gz
    """
}