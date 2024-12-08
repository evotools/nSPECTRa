// SDM processes

process sdm {
    tag "sdm"
    publishDir "${params.outdir}/sdm/raw", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/sdm-environment.yml" : null}
    

    input:
    tuple val(samplename), path(samplelist), val(chrom), path(vcf), path(tbi)
    path ancfasta
    path ancfai

    output:
    tuple val(samplename), path("sdm.${samplename}.${chrom}.txt.gz")

    script:
    """
    sdm ${vcf} ${samplelist} ${chrom} ${ancfasta} sdm.${samplename}.${chrom}
    gzip sdm.${samplename}.${chrom}.txt
    """

    stub:
    """
    touch sdm.${samplename}.${chrom}.txt.gz
    """
}

process filter_sdm {
    tag "filter_sdm"
    label "renv_large"
    publishDir "${params.outdir}/sdm/filtered", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}

    input:
    tuple val(samplename), path(samplelist), path("sdms/*")

    output:
    tuple val(samplename), path(samplelist), path("sdm.${samplename}.filtered.RData"), emit: 'rdata'
    tuple val(samplename), path(samplelist), path("sdm.${samplename}.filtered.bed"), emit: 'bed'

    script:
    """
    sdmFilter ${samplename}
    """

    stub:
    """
    touch sdm.${samplename}.filtered.RData
    """
}

// Split data in/out repetitive elements.
process repeat_mask_split_sdm {
    tag "split_sdm"
    label "medium"
    publishDir "${params.outdir}/sdm/repeat", mode: "${params.publish_dir_mode}", overwrite: true
    

    input:
    tuple val(samplename), path(samplelist), path(inbed), path("rm.bed")

    output:
    tuple val(samplename), path(samplelist), path("${inbed.baseName}.in_repeat.bed")
    tuple val(samplename), path(samplelist), path("${inbed.baseName}.out_repeat.bed")

    script:
    """
    bedtools intersect -a ${inbed} -b rm.bed -u > ${inbed.baseName}.in_repeat.bed
    bedtools intersect -a ${inbed} -b rm.bed -v > ${inbed.baseName}.out_repeat.bed
    """

    stub:
    """
    touch ${inbed.baseName}.in_repeat.bed
    touch ${inbed.baseName}.out_repeat.bed
    """
}


process count_sdm {
    tag "count_sdm"
    label "renv_large"
    publishDir "${params.outdir}/sdm/counts/${samplelist.simpleName}", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    

    input:
    tuple val(samplename), path(samplelist), path(filtered)


    output:
    path "*.txt"
    path "indChanges.${samplename}.RData"

    script:
    """
    echo "Run counter"
    getCounts3 ${samplename} ${samplelist} ${filtered}
    echo
    echo "Run individuals"
    getIndchange ${samplename} ${samplelist} 
    """

    stub:
    """
    touch indChanges.${samplename}.RData
    touch changeCounts.${samplename}.txt
    touch changeCounts_uniqued.${samplename}.txt
    touch doubleCounts.${samplename}.txt
    touch doubleCounts_uniqued.${samplename}.txt
    """
}


process make_ksfs {
    tag "ksfs"
    label "renv_large"
    publishDir "${params.outdir}/sdm/ksfs/${breedname}", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    

    input:
    tuple val(breedname), path(breedfile), path("sdms/*")


    output:
    path "sdm_${breedname}.txt"
    path "mnp_${breedname}.txt"
    path "disnps_${breedname}.txt"

    script:
    """
    SDMtoKsfs.R ${breedname}
    """

    stub:
    """
    touch sdm_${breedname}.txt
    touch mnp_${breedname}.txt
    touch disnps_${breedname}.txt
    """
}


process sdm_plot {
    tag "sdm_plot"
    label "renv_large"
    publishDir "${params.outdir}/sdm/plot/", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    

    input:
    //tuple val(breedname), path(breedfile)
    path raw_sdms


    output:
    path "plot_sdm_mutSpectra.pdf"

    script:
    """
    sdmCounts.R ./
    """

    stub:
    """
    touch plot_sdm_mutSpectra.pdf
    """
}

process sdm_matrix {
    tag "sdm_plot"
    label "renv_large"
    publishDir "${params.outdir}/sdm/matrix/", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    
    input:
    path inputs

    output:
    path "sdm_${params.species.capitalize()}_K3.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    ## Import SDMs
    allfiles <- list.files('./', pattern='*.txt', recursive = T, full.names = T)
    singlecounts <- allfiles[grepl('doubleCounts\\\\.', allfiles)]

    ## Load and collate all data
    myfiles<-list()
    for(f in singlecounts) {
        fname = str_match(f, "doubleCounts\\\\.*(.*?)\\\\.txt")[,2]
        print(fname)
        myfiles[[fname]]<-read_tsv(f, col_names = TRUE) %>% mutate(Species = '${params.species.capitalize()}')
    }

    ## Combine files
    dat<-bind_rows(myfiles, .id='Breed')
    dat <- dat[nchar(dat\$Codon1) == 3 & nchar(dat\$Codon2) == 3 & nchar(dat\$Codon3) == 3, ]
    dat <- dat[!grepl('N', dat\$Codon1) & !grepl('N', dat\$Codon2) & !grepl('N', dat\$Codon3), ]

    ## Generate counts
    counts = dat %>%
        select(Species, Breed, Ind, Codon1, Codon2, Codon3, Count) %>%
        unite('IID', c('Species', 'Breed', 'Ind'), sep = '-') %>%
        unite('Change', c('Codon1', 'Codon2', 'Codon3'), sep = '>') %>%
        group_by(IID, Change) %>%
        summarise(Count = sum(Count)) %>%
        pivot_wider(values_from = Count, names_from = Change)
    counts[is.na(counts)] = 0
    counts<-counts %>%
        separate(IID, c('Species', "Breed", "sample"), extra = "merge", sep = '-') %>%
        select(-c("Species", "Breed"))
    write.table(counts, file = "sdm_${params.species.capitalize()}_K3.csv", sep = ";", row.names = F, col.names = T, quote=F)
    """

    stub:
    """
    touch sdm_${params.species.capitalize()}_K3.csv
    """
}