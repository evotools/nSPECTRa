
/*
 * Phase 7: run mutyper
 */

process mutyper_variant {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/chrom_res", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(chrom), val(start), val(end), path(vcf), path(tbi), val(k)
    path ancfasta
    path ancfai
    path masks_ch

    output:
    tuple val(k), 
        val(chrom), val(start), val(end), 
        path("mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz"), 
        path("mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz.tbi")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz
    touch mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz.tbi
    """

    script:
    if (params.annotation)
    """
    echo "Run mutyper (variants)"
    bcftools view --threads ${task.cpus} -v snps -r ${chrom}:${start}-${end} -m2 -M2 ${vcf} | \
        bedtools intersect -header -v -a - -b ${masks_ch} | \
        sed 's/_pilon//g' | \
        vcffixup - | \
        mutyper variants --k ${k} --strand_file ${params.annotation} ${ancfasta} - | \
        bgzip -c > mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz &&
        tabix -p vcf mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz
    """
    else
    """
    echo "Run mutyper (variants)"
    bcftools view --threads ${task.cpus} -v snps -r ${chrom}:${start}-${end} -m2 -M2 ${vcf} | \
        bedtools intersect -header -v -a - -b ${masks_ch} | \
        sed 's/_pilon//g' | \
        vcffixup - | \
        mutyper variants --k ${k} ${ancfasta} - | \
        bgzip -c > mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz &&
        tabix -p vcf mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.vcf.gz
    """
}

process mutyper_spectra {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/chrom_res", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(k),
        val(chrom), val(start), val(end), 
        path(vcf),
        path(tbi)

    output:
    tuple val(k),
        val(chrom), val(start), val(end),
        path("mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.txt")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.txt
    """

    script:
    """
    echo "Run mutyper (variants)"
    mutyper spectra ${vcf} > mutationSpectra_${params.reference}_${chrom}_${start}-${end}_${k}.txt
    """
}

process mutyper_concat {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/chrom_res", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    tuple val(k),
        val(chrom), val(start), val(end), 
        path("vcfs/*"),
        path("vcfs/*")

    output:
    tuple val(k),
        path("mutyper_${params.reference}_${k}.vcf.gz"),
        path("mutyper_${params.reference}_${k}.vcf.gz.tbi")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${k}.txt
    """

    script:
    """
    echo "Run mutyper (variants)"
    bcftools concat -O u vcfs/*.vcf.gz | \
        bcftools sort -O z > mutyper_${params.reference}_${k}.vcf.gz
    bcftools index -t mutyper_${params.reference}_${k}.vcf.gz
    """
}


process count_mutations {
    tag "medium"
    label "medium"

    input:
    tuple val(k), val(chrom), val(start), val(end), path(vcf), path(tbi), path(levels)

    output:
    tuple val(k), path("mutationSpectra_${params.reference}_${k}_*.tsv")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${k}_${chrom}_${start}-${end}.tsv
    """

    script:
    """
    compute_spectra -i ${vcf} -k ${levels} -o mutationSpectra_${params.reference}_${k}_${chrom}_${start}-${end}.tsv
    """
}


process combine_counts {
    tag "count_mutations"
    label "medium"
    cpus 1
    publishDir "${params.outdir}/mutyper/full_counts", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path("tsvs/*")

    output:
    tuple val(k), path("mutationSpectra_${params.reference}_${k}.tsv")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${k}.tsv
    """

    script:
    """
    combine_matrix -i ./tsvs/ -o mutationSpectra_${params.reference}_${k}.tsv
    """
}

process count_mutations_csq {
    tag "count_mutations"
    label "medium_largemem"

    input:
    tuple val(k), val(chrom), val(start), val(end), path(vcf), path(tbi), path(levels), path(priority)

    output:
    tuple val(k), path("mutationSpectra_${params.reference}_${k}_*.csq.tsv")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${k}_${chrom}_${start}-${end}.csq.tsv
    """

    script:
    """
    compute_spectra_class -i ${vcf} -k ${levels} -c ${priority} -o mutationSpectra_${params.reference}_${k}_${chrom}_${start}-${end}.csq.tsv
    """
}


process combine_csqs {
    tag "count_mutations"
    label "medium_largemem"
    cpus 1
    publishDir "${params.outdir}/mutyper/full_counts_csq", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path("tsvs/*")

    output:
    tuple val(k), path("mutationSpectra_${params.reference}_${k}.csq.tsv")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${k}.csq.tsv
    """

    script:
    """
    combine_matrix -i ./tsvs/ -o mutationSpectra_${params.reference}_${k}.csq.tsv
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
    tuple val(k), path("mutyper_mutationSpectra_${params.reference}_${k}.csv")
    
    stub:
    """
    touch mutyper_mutationSpectra_${params.reference}_${k}.csv
    """

    script:
    $/
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

    write.csv2(mutSpectra, "mutyper_mutationSpectra_${params.reference}_${k}.csv")
    /$
}


process group_results {
    tag "group_results"
    label "medium"
    publishDir "${params.outdir}/mutyper/results", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), val(chrom), val(start), val(end), path(data)

    output:
    tuple val(k), path("mutyper_mutationSpectra_${params.reference}_${k}.csv")

    
    stub:
    """
    touch mutyper_mutationSpectra_${params.reference}_${k}.csv
    """

    script:
    """
    CombineMutyper ${k} > mutyper_mutationSpectra_${params.reference}_${k}.csv
    """
}

process plot_results {
    tag "group_results"
    label "renv"
    publishDir "${params.outdir}/mutyper/plots", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(spectra)

    output:
    tuple val(k), path("plot_mutyper_mutSpectra_${params.reference}_${k}.pdf")

    
    stub:
    """
    touch plot_mutyper_mutSpectra_${params.reference}_${k}.pdf
    """

    script:
    $/
    #!/usr/bin/env Rscript
    options(stringsAsFactors = F, warn=-1, message = FALSE, readr.num_columns = 0, dplyr.summarise.inform = FALSE)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(reshape2, quietly = TRUE))
    suppressPackageStartupMessages(library(ggfortify, quietly = TRUE))
    suppressPackageStartupMessages(library(ggforce, quietly = TRUE))
 
    mutSpectra = read_csv2("${spectra}") 

    mutSpectra[,-1]<-mutSpectra[,-1]/rowSums(mutSpectra[,-1], na.rm=T)
    mutSpectra<-mutSpectra %>% separate(sample, c("Breed", "Id"), extra = "merge")
    pca_res <- prcomp(mutSpectra[,-c(1,2)], scale. = TRUE)
    pdf("plot_mutyper_mutSpectra_${params.reference}_${k}.pdf", height = 16, width = 16)
    autoplot(pca_res, data=mutSpectra, colour = 'Breed') + geom_mark_ellipse()
    dev.off()
    /$
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

    
    stub:
    """
    touch ksfs_${samplename}_${k}.tsv
    """

    script:
    """
    echo "Run mutyper (ksfs)"
    bcftools view --threads ${task.cpus} -S ${samplelist} ${vcf} | \
        mutyper ksfs - > ksfs_${samplename}_${k}.tsv
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
    path "K${k}_counts.txt"

    stub:
    """
    touch K${k}_counts.txt
    """

    script:
    """
    jellyfish count -t 4 -m ${k} -s 3G ${ancfa} -o K${k}.jf
    jellyfish dump -o K${k}_counts.txt -c K${k}.jf && rm K${k}.jf
    """
}

process normalize_results {
    tag 'kcnt'
    label 'medium'
    publishDir "${params.outdir}/mutyper/results_corrected", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(k), path(counts)
    path normalizers

    output:
    path "${counts.baseName}.*.csv"
    path "${counts.baseName}.*.csv"

    stub:
    """
    touch ${counts.baseName}.Knorm.csv
    touch ${counts.baseName}.KCnorm.csv
    """

    script:
    """
    CORRECT_COUNTS -s ${counts} -k K${k}_counts.txt -m 1 
    CORRECT_COUNTS -s ${counts} -k K${k}_counts.txt -m 2 
    """
}