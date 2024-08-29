process sdm {
    tag "sdm"
    publishDir "${params.outdir}/sdm/raw", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/sdm-environment.yml" : null}
    

    input:
    path vcf
    path tbi
    path reffasta
    path reffai
    tuple val(samplename), path(samplelist), val(idx), val(contig)

    output:
    path "sdm.${samplename}.${contig}.txt.gz"

    stub:
    """
    touch sdm.${samplename}.${contig}.txt.gz
    """

    script:
    """
    sdm ${vcf} ${samplelist} ${contig} ${reffasta} sdm.${samplename}.${contig}
    gzip sdm.${samplename}.${contig}.txt
    """
}

process filter_sdm {
    tag "filter_sdm"
    label "renv_large"
    publishDir "${params.outdir}/sdm/filtered", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    

    input:
    tuple val(samplename), path(samplelist)
    path sdms

    output:
    tuple val(samplename), path(samplelist), path("sdm.${samplename}.filtered.RData"), emit: 'rdata'
    tuple val(samplename), path(samplelist), path("sdm.${samplename}.filtered.bed"), emit: 'bed'

    stub:
    """
    touch sdm.${samplename}.filtered.RData
    """

    script:
    """
    sdmFilter ${samplename}
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

    stub:
    """
    touch ${inbed.baseName}.in_repeat.bed
    touch ${inbed.baseName}.out_repeat.bed
    """

    script:
    """
    bedtools intersect -a ${inbed} -b rm.bed -u > ${inbed.baseName}.in_repeat.bed
    bedtools intersect -a ${inbed} -b rm.bed -v > ${inbed.baseName}.out_repeat.bed
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

    stub:
    """
    touch indChanges.${samplename}.RData
    touch changeCounts.${samplename}.txt
    touch changeCounts_uniqued.${samplename}.txt
    touch doubleCounts.${samplename}.txt
    touch doubleCounts_uniqued.${samplename}.txt
    """

    script:
    """
    echo "Run counter"
    getCounts3 ${samplename} ${samplelist} ${filtered}
    echo
    echo "Run individuals"
    getIndchange ${samplename} ${samplelist} 
    """
}


process make_ksfs {
    tag "ksfs"
    label "renv_large"
    publishDir "${params.outdir}/sdm/ksfs/${breedname}", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null}
    

    input:
    tuple val(breedname), path(breedfile)
    path raw_sdms


    output:
    path "sdm_${breedname}.txt"
    path "mnp_${breedname}.txt"
    path "disnps_${breedname}.txt"

    stub:
    """
    touch sdm_${breedname}.txt
    touch mnp_${breedname}.txt
    touch disnps_${breedname}.txt
    """

    script:
    """
    SDMtoKsfs.R ${breedname}
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

    stub:
    """
    touch plot_sdm_mutSpectra.pdf
    """

    script:
    """
    sdmCounts.R ./
    """
}