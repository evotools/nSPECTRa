
/*
 * Phase 7: run mutyper
 */

process mutyper {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/chrom_res", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    path vcf
    path tbi
    path ancfasta
    path ancfai
    path masks_ch
    tuple val(idx), val(region), val(k)

    output:
    tuple val(k), val(region), path("mutationSpectra_${params.reference}_${region}_${k}.txt")

    
    stub:
    """
    touch mutationSpectra_${params.reference}_${region}_${k}.txt
    """

    script:
    if (params.annotation)
    """
    echo "Run mutyper (variants)"
    bcftools view -v snps  -r ${region} -m2 -M2 ${vcf} |
        bedtools intersect -header -v -a - -b ${masks_ch} |
        sed 's/_pilon//g' |
        vcffixup - | 
        mutyper variants --k ${k} --strand_file ${params.annotation} ${ancfasta} - |
        mutyper spectra - > mutationSpectra_${params.reference}_${region}_${k}.txt
    """
    else
    """
    echo "Run mutyper (variants)"
    bcftools view -v snps  -r ${region} -m2 -M2 ${vcf} |
        bedtools intersect -header -v -a - -b ${masks_ch} | \
        sed 's/_pilon//g' |
        vcffixup - | 
        mutyper variants --k ${k} ${ancfasta} - |
        mutyper spectra - > mutationSpectra_${params.reference}_${region}_${k}.txt
    """
}


process mutyper_full {
    tag "mutyper_full"
    label "medium"
    publishDir "${params.outdir}/mutyper/vcfs", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    path vcf
    path tbi
    path ancfasta
    path ancfai
    path masks_ch
    val k

    output:
    tuple val(k), path("mutyper_${params.reference}_${k}.vcf.gz"), path("mutyper_${params.reference}_${k}.vcf.gz.tbi")
    
    stub:
    """
    touch mutyper_${params.reference}_${k}.vcf.gz
    touch mutyper_${params.reference}_${k}.vcf.gz.tbi
    """

    script:
    if (params.annotation)
    """
    echo "Run mutyper (variants)"
    bcftools view -v snps -m2 -M2 ${vcf} |
        bedtools intersect -header -v -a - -b ${masks_ch} | 
        sed 's/_pilon//g' |
        vcffixup - | 
        mutyper variants --k ${k} --strand_file ${params.annotation} ${ancfasta} - |
        bgzip -c > mutyper_${params.reference}_${k}.vcf.gz
    tabix -p vcf mutyper_${params.reference}_${k}.vcf.gz
    """
    else
    """
    echo "Run mutyper (variants)"
    bcftools view -v snps -m2 -M2 ${vcf} |
        bedtools intersect -header -v -a - -b ${masks_ch} |  
        sed 's/_pilon//g' |
        vcffixup - | 
        mutyper variants --k ${k} ${ancfasta} - |
        bgzip -c > mutyper_${params.reference}_${k}.vcf.gz
    tabix -p vcf mutyper_${params.reference}_${k}.vcf.gz
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
    tuple val(k), val(region), path(data)

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
 
    mutSpectra = read_csv2("${spectra}") 

    mutSpectra[,-1]<-mutSpectra[,-1]/rowSums(mutSpectra[,-1], na.rm=T)
    mutSpectra<-mutSpectra %>% separate(sample, c("Breed", "Id"), extra = "merge")
    pca_res <- prcomp(mutSpectra[,-c(1,2)], scale. = TRUE)
    pdf("plot_mutyper_mutSpectra_${params.reference}_${k}.pdf", height = 8, width = 12)
    autoplot(pca_res, data=mutSpectra, colour = 'Breed')
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
    bcftools view -S ${samplelist} ${vcf} |
        mutyper ksfs - > ksfs_${samplename}_${k}.tsv
    """
}
