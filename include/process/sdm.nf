process sdm {
    tag "sdm"
    publishDir "${params.outdir}/sdm/raw", mode: "${params.publish_dir_mode}", overwrite: true
    

    input:
    file vcf
    file tbi
    file reffasta
    file reffai
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
    

    input:
    tuple val(samplename), file(samplelist)
    path sdms

    output:
    tuple val(samplename), path(samplelist), path("sdm.${samplename}.filtered.RData")

    stub:
    """
    touch sdm.${samplename}.filtered.RData
    """

    script:
    """
    sdmFilter ${samplename}
    """
}


process count_sdm {
    tag "count_sdm"
    label "renv_large"
    publishDir "${params.outdir}/sdm/counts/${samplelist.simpleName}", mode: "${params.publish_dir_mode}", overwrite: true
    

    input:
        tuple val(samplename), path(samplelist), path(filtered)


    output:
    path "*.txt"
    path "indChanges.*"

    stub:
    """
    touch indChanges.${samplename}.Rdata
    touch changeCounts.${samplename}.txt
    touch changeCounts_uniqued.${samplename}.txt
    touch doubleCounts.${samplename}.txt
    touch doubleCounts_uniqued.${samplename}.txt
    """

    script:
    """
    echo "Run counter"
    getCounts3 ${samplename} ${samplelist} ${filtered}
    echo "Run individuals"
    getIndchange ${samplename} ${samplelist} 
    """
}

