


/*
 * Compute constrained elements   
 */

process hal4d {
    tag "h4d"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "large_onecore"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path hal
    path exon_bed

    output:
    tuple path(hal), path("neutralRegions.bed")
    
    
    stub:
    """
    touch neutralRegions.bed 
    """

    script:
    """
    hal4dExtract --hdf5InMemory ${hal} "${params.reference}" ${exon_bed} neutralRegions.bed 
    """
}

process halTree {
    tag "pFit"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path hal

    output:
    env TREE
    
    
    stub:
    """
    TREE="((spp1, spp2), spp3)"
    """

    script:
    """
    TREE=`halStats --tree ${hal}`
    """
}

process phyloFit {
    tag "pFit"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    tuple path(ss), val(TREE)

    output:
    path "neutralModel.mod"
    
    
    stub:
    """
    touch neutralModel.mod
    """

    script:
    """
    phyloFit --tree "${TREE}" --msa-format SS --subst-mod SSREV --sym-freqs --precision HIGH --out-root neutralModel ${ss} 
    """
}


process make4dmaf {
    tag "maf"
    publishDir "${params.outdir}/PHYLOP/MAF", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_vlargemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    tuple path(hal), path(bedfile), val(GENOMES)

    output:
    tuple path(hal), path("4d.maf"), val(GENOMES)
    
    stub:
    """
    touch 4d.maf
    """

    script:
    """
    hal2mafMP.py \
        ${hal} \
        4d.maf \
        --noDupes \
        --targetGenomes "${GENOMES}" \
        --numProc ${task.cpus} \
        --refGenome ${params.reference} \
        --refTargets ${bedfile} \
        --hdf5InMemory
    sed -i -e 2d 4d.maf
    """
}

process hal_genomes {
    tag "msa"
    publishDir "${params.outdir}/PHYLOP/MSA", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path hal

    output:
    env(GENOMES)
    
    
    stub:
    """
    GENOMES='spp1,spp2,spp3'
    """

    script:
    """
    GENOMES=\$( halStats ${hal} --genomes | sed 's/ /\\n/g' | paste -sd, )
    """
}

process msa_view {
    tag "msa"
    publishDir "${params.outdir}/PHYLOP/MSA", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    tuple path(hal), path(maf), val(GENOMES)

    output:
    path "${maf.simpleName}.ss"
    
    
    stub:
    """
    touch ${maf.simpleName}.ss
    """

    script:
    """
    MAF_RENAME -m ${maf} -i "${GENOMES}" -o renamed.maf
    CONV=\$( cat conversion.csv )
    msa_view -o SS -z --in-format MAF --aggregate "\$CONV" renamed.maf | \
        sed "s/NAMES = \$CONV/NAMES = ${GENOMES}/g" > ${maf.simpleName}.ss
    """
}


process phyloPtrain {
    tag "pptrain"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_multi"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path hal
    path ss

    output:
    path "neutralModel.mod"
    
    
    stub:
    """
    touch neutralModel.mod
    """

    script:
    """
    halPhyloPTrain.py ${hal} ${params.reference} ${neutral} neutralModel.mod --numProc ${task.cpus} 
    """
}


process phyloP {
    tag "ppmp"
    publishDir "${params.outdir}/PHYLOP/WIG/", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    tuple val(n), val(chr), path(hal), path(model)

    output:
    tuple val(n), val(chr), path(hal), path(model), path("phylop_${chr}.wig"), path("${params.reference}.sizes")
    
    stub:
    """
    touch phylop_${chr}.wig
    """

    script:
    """
    halPhyloPMP.py \
        ${hal} \
        ${params.reference} \
        ${model} \
        phylop_${chr}.wig \
        --refSequence ${chr} \
        --chromSizes ${params.reference}.sizes \
        --numProc ${task.cpus} 
    """
}


process wig2bedgraph {
    tag "ppmp"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    afterScript "rm ${wig.baseName}.bw"

    input:
    tuple val(n), val(chr), path(hal), path(model), path(wig), path(sizes)

    output:
    path "${wig.baseName}.bed"
    
    
    stub:
    """
    touch ${wig.baseName}.bed
    """

    script:
    """
    wigToBigWig ${wig} ${sizes} ${wig.baseName}.bw
    bigWigToBedGraph ${wig.baseName}.bw  /dev/stdout | \
        sort -k1,1 -k2,2n --parallel ${task.cpus} - > ${wig.baseName}.bed
    """
}


process bedtobigwig {
    tag "bw2bed"
    publishDir "${params.outdir}/PHYLOP/BED/CHR", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda { params.enable_conda ? "ucsc-bedgraphtobigwig" : null }
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path bw

    output:
    path "${bw.simpleName}.bed"
    
    
    stub:
    """
    touch ${bw.simpleName}.bed
    """

    script:
    """
    bigWigToBedGraph ${bw} /dev/stdout | bedtools sort -i - > ${bw.simpleName}.bed
    """
}

process combine_bed {
    tag "bed"
    publishDir "${params.outdir}/PHYLOP/BED", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    path beds

    output:
    path "phylop.bed"
    
    
    stub:
    """
    touch phylop.bed
    """

    script:
    """
    cat ${beds} | bedtools sort -i - > phylop.bed
    """
}


process extract_conserved {
    publishDir "${params.outdir}/PHYLOP/CONSERVED", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    path "phylop.bed"

    output:
    path "phylop_conserved.bed"
    
    
    stub:
    """
    touch phylop_conserved.bed
    """

    script:
    """
    #!/usr/bin/env Rscript
    options(scipen = 999)
    library(tidyverse)

    # Fetch input BEDgraph
    bed <- read_table('phylop.bed', col_names = c('chrom', 'start', 'end', 'phylop'))

    # Define Q3 and IQR
    q3 <- as.vector(quantile(bed\$phylop, probs=c(0.75)))
    iqr <- IQR(bed\$phylop)

    # Extract conserved elements as positive outliers
    conserved <- bed %>% filter( phylop > q3+iqr*1.5 )

    # Save conserved regions
    write.table(
        conserved,
        'phylop_conserved.bed',
        sep = '\t',
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
    )
    """
}


process vcf_drop_conserved {
    tag "filt"
    publishDir "${params.outdir}/PHYLOP/VCF", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium"

    input:
    path vcf
    path tbi
    path bed

    output:
    path "${vcf.simpleName}.non-conserved.vcf.gz", emit: vcf
    path "${vcf.simpleName}.non-conserved.vcf.gz.tbi", emit: tbi
    
    
    stub:
    """
    touch ${vcf.simpleName}.non-conserved.vcf.gz
    touch ${vcf.simpleName}.non-conserved.vcf.gz.tbi
    """

    script:
    """
    bedtools intersect -header -v -a ${vcf} -b ${bed} | bgzip -c > ${vcf.simpleName}.non-conserved.vcf.gz
    tabix -p vcf ${vcf.simpleName}.non-conserved.vcf.gz
    """
}

