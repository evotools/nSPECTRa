
/*
 * Phase 1. From hal to ancestral fasta  
 */
process hal2maf {
    tag "hal2maf"
    publishDir "${params.outdir}/MAF", mode: "${params.publish_dir_mode}", overwrite: true

    label "largemem"

    input:
    path HAL
    path CACTUS

    output:
    path "${params.reference}_${params.target}.maf"  
    
    
    stub:
    """
    touch ${params.reference}_${params.target}.maf
    """

    script:
    """
    ${CACTUS}/bin/hal2maf \
        --refGenome ${params.reference} \
        --targetGenomes ${params.target} \
        --hdf5InMemory ${HAL} ${params.reference}_${params.target}.maf
    """
}


process hal2chain {
    tag "hal2chain"
    label "largemem"
    publishDir "${params.outdir}/CHAIN", mode: "${params.publish_dir_mode}", overwrite: true

    
    input:
    path reference
    path target
    path HAL
    path CACTUS

    output:
    path "${params.target}_${params.reference}.chain"  
    
    
    stub:
    """
    touch ${params.target}_${params.reference}.chain
    """

    script:
    """
    ${CACTUS}/bin/hal2maf \
        --refGenome ${params.target} \
        --targetGenomes ${params.reference} \
        --hdf5InMemory ${HAL} stdout | \
            maf-convert chain - > ${params.target}_${params.reference}.chain
    """
}


process halSnps {
    tag "halsnps"
    label "largemem"
    publishDir "${params.outdir}/SNPs", mode: "${params.publish_dir_mode}", overwrite: true

    
    input:
    path HAL
    path CACTUS

    output:
    path "SNPs_${params.reference}_${params.target}.tsv"  
    
    
    stub:
    """
    touch SNPs_${params.reference}_${params.target}.tsv 
    """

    script:
    """
    ${CACTUS}/bin/halSnps ${HAL} ${params.reference} ${params.target} \
        --tsv SNPs_${params.reference}_${params.target}.tsv --hdf5InMemory && \
        rm ./\${halname}.hal
    """
}