


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

    script:
    """
    hal4dExtract --hdf5InMemory ${hal} "${params.reference}" ${exon_bed} neutralRegions.bed 
    """
    
    stub:
    """
    touch neutralRegions.bed 
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
    env "TREE"

    script:
    """
    TREE=`halStats --tree ${hal}`
    """

    stub:
    """
    TREE="((spp1, spp2), spp3)"
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

    script:
    """
    phyloFit --tree "${TREE}" --msa-format SS --subst-mod SSREV --sym-freqs --precision HIGH --out-root neutralModel ${ss} 
    """
    
    stub:
    """
    touch neutralModel.mod
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

    script:
    if (params.cactus_hal2maf)
    """
    export HOME=\$PWD
    toil config test.yaml
    mkdir toil-work
    mkdir toil-coord
    cactus-hal2maf \
        ./js \
        ${hal} \
        4d.maf \
        --config test.yaml \
        --refGenome ${params.reference} \
        --bedRanges ${bedfile} \
        --targetGenomes "${GENOMES}" \
        --batchCores ${task.cpus} \
        --defaultCores ${task.cpus} \
        --workDir ./toil-work \
        --coordinationDir ./toil-coord \
        --cleanWorkDir onSuccess
    """
    else
    """
    hal2maf \
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
    
    stub:
    """
    touch 4d.maf
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
    env "GENOMES"

    script:
    """
    GENOMES=\$( halStats ${hal} --genomes | sed 's/ /\\n/g' | paste -sd, )
    """
    
    stub:
    """
    GENOMES='spp1,spp2,spp3'
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

    script:
    """
    MAF_RENAME -m ${maf} -i "${GENOMES}" -o renamed.maf
    CONV=\$( cat conversion.csv )
    msa_view -o SS -z --in-format MAF --aggregate "\$CONV" renamed.maf | \
        sed "s/NAMES = \$CONV/NAMES = ${GENOMES}/g" > ${maf.simpleName}.ss
    """
    
    stub:
    """
    touch ${maf.simpleName}.ss
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

    script:
    """
    halPhyloPTrain.py ${hal} ${params.reference} ${ss} neutralModel.mod --numProc ${task.cpus} 
    """
    
    stub:
    """
    touch neutralModel.mod
    """
}


process phyloP {
    tag "ppmp"
    publishDir "${params.outdir}/PHYLOP/WIG/", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    tuple path(hal), path(BED)
    path model

    output:
    tuple path("phylop_${chr}.wig"), path("${params.reference}.sizes")

    script:
    """
    halPhyloPMP.py \
        ${hal} \
        ${params.reference} \
        ${model} \
        ${BED.simpleName}.wig \
        --refBed ${BED} \
        --chromSizes ${params.reference}.sizes \
        --numProc ${task.cpus} 
    """
    
    stub:
    """
    touch phylop_${chr}.wig
    """
}


process wig2bedgraph {
    tag "ppmp"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}
    afterScript "rm ${wig.baseName}.bw"

    input:
    tuple path(wig), path(sizes)

    output:
    path "${wig.baseName}.bed"

    script:
    """
    wigToBigWig ${wig} ${sizes} ${wig.baseName}.bw
    bigWigToBedGraph ${wig.baseName}.bw  /dev/stdout | \
        sort -k1,1 -k2,2n --parallel ${task.cpus} - > ${wig.baseName}.bed
    """
    
    stub:
    """
    touch ${wig.baseName}.bed
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

    script:
    """
    bigWigToBedGraph ${bw} /dev/stdout | bedtools sort -i - > ${bw.simpleName}.bed
    """
    
    stub:
    """
    touch ${bw.simpleName}.bed
    """
}

process combine_bed {
    tag "bed"
    publishDir { "${params.outdir}/${outdir}/BED" }, mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    path beds
    val outname
    val outdir

    output:
    path "${outname}.bed"

    script:
    """
    cat ${beds} | sort -k 1,1 -k2,2n --parallel ${task.cpus} - > ${outname}.bed
    """
    
    stub:
    """
    touch phylop.bed
    """
}


process extract_conserved {
    publishDir "${params.outdir}/PHYLOP/CONSERVED", mode: "${params.publish_dir_mode}", overwrite: true
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    path "phylop.bed"

    output:
    path "phylop_conserved.bed"

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
    
    stub:
    """
    touch phylop_conserved.bed
    """
}


process vcf_drop_intervals {
    tag "filt"
    publishDir "${params.outdir}/PHYLOP/VCF", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium"

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path bed
    val tag

    output:
    path "${vcf.simpleName}.${tag}.vcf.gz", emit: vcf
    path "${vcf.simpleName}.${tag}.vcf.gz.tbi", emit: tbi

    script:
    """
    bedtools intersect -header -v -a ${vcf} -b ${bed} | \
        bgzip -@ ${task.cpus > 1 ? task.cpus - 1 : 1} -c > ${vcf.simpleName}.${tag}.vcf.gz
    tabix -p vcf ${vcf.simpleName}.${tag}.vcf.gz
    """
    
    stub:
    """
    touch ${vcf.simpleName}.${tag}.vcf.gz
    touch ${vcf.simpleName}.${tag}.vcf.gz.tbi
    """
}

// Processes for the phastBias methods
process GENOME_INTERVALS {
    conda "bioconda::pysam=0.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py39hcada746_0' :
        'quay.io/biocontainers/pysam:0.22.1--py39hcada746_0' }"

    input:
    path sizes

    output:
    path "intervals_*.bed"

    script:
    """
    #!/usr/bin/env python
    import pysam
    import re

    sizes = open("${sizes}")
    n = 1
    target_size = ${params.chunk_size}
    proc_size = 0
    tmp_list = []
    for line in sizes:
        seq_id, seq_len = line.strip().split()
        seq_len = int(seq_len)
        if seq_len > target_size:
            for i in range(0, seq_len, target_size):
                with open(f"intervals_{n}.bed", "w") as bedfile:
                    start = i
                    end = min(i + target_size, seq_len)
                    bedfile.write(f'{seq_id}\\t{start}\\t{end}\\n')
                    n+=1
        else:
            with open(f"intervals_{n}.bed", "w") as bedfile:
                bedfile.write(f'{seq_id}\\t0\\t{seq_len}\\n')
                n+=1
    if len(tmp_list) > 0:
        for line in tmp_list:
            bedfile.write(line)
    """
}

process create_maf {
    tag "hal2maf"
    publishDir "${params.outdir}/MAF", mode: "${params.publish_dir_mode}", overwrite: true
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    label "largemem"

    input:
    tuple path(HAL), path(BED)

    output:
    tuple path("${BED.simpleName}.maf"), path(BED)
    

    script:
    if (params.cactus_hal2maf)
    """
    export HOME=\$PWD
    mkdir toil-work
    mkdir toil-coord
    toil config test.yaml
    cactus-hal2maf \
        ./js \
        ${HAL} \
        ${BED.simpleName}.maf \
        --config test.yaml \
        --refGenome ${params.reference} \
        --noAncestors \
        --bedRanges ${BED} \
        --batchCores 1 \
        --outType single \
        --defaultCores 1 \
        --workDir ./toil-work \
        --coordinationDir ./toil-coord \
        --cleanWorkDir onSuccess
    """
    else
    """
    hal2maf \
        --noAncestors \
        --refGenome ${params.reference} \
        --targetGenomes "${GENOMES}" \
        --refTargets ${BED} \
        --hdf5InMemory ${HAL} alignments.${CHROM}.maf
    """
    
    stub:
    """
    touch alignments.${CHROM}.maf
    """
}



process phastBias {
    tag "phastBias"
    publishDir "${params.outdir}/PHAST/gBGC", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda {params.enable_conda ? "${baseDir}/envs/phast_environment.yml" : null}

    input:
    tuple path(maf), path(BED)
    path model

    output:
    path "${maf.simpleName}.wig", emit: wig
    path "${maf.simpleName}.tracts.bed", emit: tracts_gff
    path "${maf.simpleName}.informative.bed", emit: informative_gff
    path "${maf.simpleName}.tracts.bed", emit: tracts_bed
    path "${maf.simpleName}.informative.bed", emit: informative_bed

    script:
    """
    CHROM=\$( head -1 $BED | cut -f1 )
    phastBias \
        --informative-fn ${maf.simpleName}.informative.gff \
        --output-tracts ${maf.simpleName}.tracts.gff \
        ${maf} \
        ${model} \
        ${params.reference} | \
        sed "s/chrom=${params.reference}/chrom=\$CHROM/" > ${maf.simpleName}.wig
    awk -v var=\$CHROM '{OFS="\\t"; print var, \$4, \$5, \$3}' ${maf.simpleName}.informative.bed > ${maf.simpleName}.informative.bed
    awk -v var=\$CHROM '{OFS="\\t"; print var, \$4, \$5, \$3}' ${maf.simpleName}.tracts.bed > ${maf.simpleName}.tracts.bed
    """
    
    stub:
    """
    touch ${maf.simpleName}.wig
    touch ${maf.simpleName}.gff
    """
}

process halSize {
    tag "medium_mem"
    publishDir "${params.outdir}/MAF", mode: "${params.publish_dir_mode}", overwrite: true
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    label "largemem"

    input:
    path HAL

    output:
    path "${params.reference}.sizes"
    

    script:
    """
    halStats --sequenceStats ${params.reference} ${HAL} | \
        awk 'BEGIN{FS=","}; NR>1 && \$1!="" {print \$1"\\t"\$2}' > ${params.reference}.sizes
    """
    
    stub:
    """
    touch ${params.reference}.sizes
    """
}

process bgcFilter {
    tag "medium_mem"
    publishDir "${params.outdir}/BGC", mode: "${params.publish_dir_mode}", overwrite: true
    label "small"

    input:
    path BED

    output:
    path "${BED.simpleName}.bgc${params.bgc_threshold}.bed"
    

    script:
    """
    awk '\$4 >= ${params.bgc_threshold}' ${BED} > ${BED.simpleName}.bgc${params.bgc_threshold}.bed
    """
    
    stub:
    """
    touch ${BED.simpleName}.bgc${params.bgc_threshold}.bed
    """
}

process catsort {
    publishDir { "${params.outdir}/${OUTDIR}" }, mode: "${params.publish_dir_mode}", overwrite: true
    label "small"

    input:
    path INFILES
    val OUTNAME
    val OUTDIR
    val SORTEXPR

    output:
    path "${OUTNAME}"
    

    script:
    """
    cat ${INFILES} | sort ${SORTEXPR} > ${OUTNAME}
    """
    
    stub:
    """
    touch ${OUTNAME}
    """
}