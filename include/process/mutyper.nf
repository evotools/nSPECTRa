
/*
 * Phase 7: run mutyper
 */

process mutyper_variant {
    tag "mutyper"
    label "medium"


    input:
    tuple val(chrom), path(vcf), path(tbi), val(k)
    path ancfasta
    path ancfai
    path masks_ch

    output:
    tuple val(k), 
        val(chrom),
        path("mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz"), 
        path("mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz.tbi")


    script:
    String vcftools_filter = ""
    if (params.mutyper_min_gq || params.mutyper_max_missing){
        def min_gq = params.mutyper_min_gq ? "--minGQ ${params.mutyper_min_gq}" : ""
        def max_miss = params.mutyper_max_missing ? "--max-missing ${params.mutyper_max_missing}" : ""
        vcftools_filter = "| vcftools --vcf - ${min_gq} ${max_miss} --recode --recode-INFO-all --stdout"
    }
    if (params.annotation)
    """
    echo "Run mutyper (variants)"
    bcftools view --threads ${task.cpus} -v snps -r ${chrom} -m2 -M2 ${vcf} | \
        bedtools intersect -header -v -a - -b ${masks_ch} | \
        sed 's/_pilon//g' | \
        vcffixup - ${vcftools_filter} |\
        mutyper variants --k ${k} --strand_file ${params.annotation} ${ancfasta} - | \
        bgzip -c > mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz &&
        tabix -p vcf mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz
    """
    else
    """
    echo "Run mutyper (variants)"
    bcftools view --threads ${task.cpus} -v snps -r ${chrom} -m2 -M2 ${vcf} | \
        bedtools intersect -header -v -a - -b ${masks_ch} | \
        sed 's/_pilon//g' | \
        vcffixup - ${vcftools_filter} |\
        mutyper variants --k ${k} ${ancfasta} - | \
        bgzip -c > mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz &&
        tabix -p vcf mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz
    """
    
    stub:
    """
    touch mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz
    touch mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.vcf.gz.tbi
    """
}

process mutyper_spectra {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/chrom_res", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(k),
        val(chrom), 
        path(vcf),
        path(tbi)

    output:
    tuple val(k),
        val(chrom),
        path("mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.txt")

    script:
    """
    echo "Run mutyper (variants)"
    mutyper spectra ${vcf} > mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.txt
    """
    
    stub:
    """
    touch mutationSpectra_${params.species.capitalize()}_${chrom}_${k}.txt
    """
}

process mutyper_concat {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/vcf", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(k),
        val(chrom), 
        path("vcfs/*"),
        path("vcfs/*")

    output:
    tuple val(k),
        path("mutyper_${params.species.capitalize()}_${k}.vcf.gz"),
        path("mutyper_${params.species.capitalize()}_${k}.vcf.gz.tbi")

    script:
    """
    echo "Run mutyper (variants)"
    bcftools concat -O u vcfs/*.vcf.gz | \
        bcftools sort -O z > mutyper_${params.species.capitalize()}_${k}.vcf.gz
    bcftools index -t mutyper_${params.species.capitalize()}_${k}.vcf.gz
    """
    
    stub:
    """
    touch mutationSpectra_${params.species.capitalize()}_${k}.txt
    """
}


process consequence_table {
    tag "medium"
    label "medium"
    publishDir "${params.outdir}/mutyper/csq_table", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(vcf), path(tbi)

    output:
    tuple val(k), path("mutyper_${params.species.capitalize()}_${k}.singled_csqs.counts.tsv.gz")


    script:
    """
    bcftools +split-vep \
            -f "%CHROM\\t%POS\\t%mutation_type\\t%Consequence\\t%Codons\\n" \
            ${vcf} | \
        CSQ_UNIFORM /dev/stdin | \
        bgzip -c > mutyper_${params.species.capitalize()}_${k}.singled_csqs.counts.tsv.gz
    """
    
    stub:
    """
    touch mutyper_${params.species.capitalize()}_${k}.singled_csqs.counts.tsv.gz
    """
}


process count_mutations {
    tag "medium"
    label "medium"
    publishDir "${params.outdir}/mutyper/full_counts", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(vcf), path(tbi), path(levels)

    output:
    tuple val(k), path("mutationSpectra_${params.species.capitalize()}_${k}.tsv")

    script:
    """
    compute_spectra -i ${vcf} -k ${levels} -o mutationSpectra_${params.species.capitalize()}_${k}.tsv
    """
    
    stub:
    """
    touch mutationSpectra_${params.species.capitalize()}_${k}.tsv
    """
}


process count_mutations_csq {
    tag "count_mutations"
    label "medium_largemem"
    publishDir "${params.outdir}/mutyper/full_counts_csq", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(vcf), path(tbi), path(levels), path(priority)

    output:
    tuple val(k), path("mutationSpectra_${params.species.capitalize()}_${k}.csq.tsv")

    script:
    """
    compute_spectra_class -i ${vcf} -k ${levels} -c ${priority} -o mutationSpectra_${params.species.capitalize()}_${k}.csq.tsv
    """
    
    stub:
    """
    touch mutationSpectra_${params.species.capitalize()}_${k}.csq.tsv
    """
}


//
// Group mutyper
//
process group_results_old {
    tag "group_results"
    label "renv_large"
    publishDir "${params.outdir}/mutyper/results", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), val(region), path(data)

    output:
    tuple val(k), path("mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv")

    script:
    """
    #!/usr/bin/env Rscript
    options(stringsAsFactors = F, warn=-1, message = FALSE, readr.num_columns = 0, dplyr.summarise.inform = FALSE)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(reshape2, quietly = TRUE))
    suppressPackageStartupMessages(library(ggfortify, quietly = TRUE))
    allmuts = list.files(pattern = "_${k}.txt", recursive = T)

    dat<-list()
    for(f in allmuts)
    {
        dat[[f]] <- read_tsv(f) %>% pivot_longer(!sample, names_to="Change", "Count")
    }
    dat<-bind_rows(dat, .id="File")
 
    mutSpectra = dat %>% 
        drop_na() %>% 
        group_by(sample, Change) %>% 
        summarise(sum=sum(value)) %>% 
        pivot_wider(names_from = Change, values_from = sum)

    write.csv2(mutSpectra, "mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv")
    """
    
    stub:
    """
    touch mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv
    """
}


process group_results {
    tag "group_results"
    label "medium"
    publishDir "${params.outdir}/mutyper/results", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), val(chrom), path(data)

    output:
    tuple val(k), path("mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv")

    script:
    """
    CombineMutyper ${k} > mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv
    """
    
    stub:
    """
    touch mutyper_mutationSpectra_${params.species.capitalize()}_${k}.csv
    """
}

process plot_results {
    tag "group_results"
    label "renv"
    publishDir "${params.outdir}/mutyper/plots", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(spectra)
    path "meta.tsv"

    output:
    tuple val(k), path("plot_mutyper_mutSpectra_${params.species.capitalize()}_${k}.pdf")

    script:
    """
    #!/usr/bin/env Rscript
    options(stringsAsFactors = F, warn=-1, message = FALSE, readr.num_columns = 0, dplyr.summarise.inform = FALSE)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(reshape2, quietly = TRUE))
    suppressPackageStartupMessages(library(ggfortify, quietly = TRUE))
    suppressPackageStartupMessages(library(ggforce, quietly = TRUE))
 
    mutSpectra <- read_csv2("${spectra}") %>% select(where(~ any(. != 0)))
    # Add metadata
    meta <- read_table("meta.tsv", col_names = c("pop", "sample"))
    mutSpectra <- merge(meta, mutSpectra, by = 'sample')

    mutSpectra[,-c(1, 2)]<-mutSpectra[,-c(1, 2)]/rowSums(mutSpectra[,-c(1, 2)], na.rm=T)
    pca_res <- prcomp(mutSpectra[,-c(1,2)], scale. = TRUE)
    pdf("plot_mutyper_mutSpectra_${params.species.capitalize()}_${k}.pdf", height = 16, width = 16)
    autoplot(pca_res, data=mutSpectra, colour = 'pop') + geom_mark_ellipse()
    dev.off()
    """
    
    stub:
    """
    touch plot_mutyper_mutSpectra_${params.species.capitalize()}_${k}.pdf
    """
}

// Extract Ksfs
process ksfs {
    tag "ksfs"
    label "medium"
    publishDir "${params.outdir}/mutyper/ksfs", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(samplename), path(samplelist), val(k), path(vcf), path(tbi)

    output:
    tuple val(k), path("ksfs_${samplename}_${k}.tsv")

    script:
    """
    echo "Run mutyper (ksfs)"
    bcftools query -l ${vcf} > all_samples.txt
    grep -f all_samples.txt ${samplelist} > keep.txt
    bcftools view --threads ${task.cpus} -S keep.txt ${vcf} | \
        mutyper ksfs - > ksfs_${samplename}_${k}.tsv
    """
    
    stub:
    """
    touch ksfs_${samplename}_${k}.tsv
    """
}

process kmercount {
    tag 'kcnt'
    label 'large'
    publishDir "${params.outdir}/mutyper/kmer_counts", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path ancfa
    path ancfai
    val k

    output:
    tuple val(k), path("${params.species.capitalize()}.K${k}.txt")

    script:
    """
    kmer_count ${ancfa} ${k} > ${params.species.capitalize()}.K${k}.txt
    """

    stub:
    """
    touch ${params.species.capitalize()}.K${k}.txt
    """
}

process normalize_results {
    tag 'kcnt'
    label 'medium'
    publishDir "${params.outdir}/mutyper/normalized", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(counts), path(normalizers)

    output:
    path "*.csv"

    script:
    """
    CORRECT_COUNTS -s ${counts} -k ${normalizers} -m 1 
    CORRECT_COUNTS -s ${counts} -k ${normalizers} -m 2 
    """

    stub:
    """
    touch ${counts.baseName}.Knorm.csv
    touch ${counts.baseName}.KCnorm.csv
    """
}