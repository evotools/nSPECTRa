


/*
 * Compute constrained elements   
 */

process hal4d {
    tag "h4d"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "large_onecore"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path HAL

    output:
    path "neutralRegions.bed"
    
    
    stub:
    """
    touch neutralRegions.bed 
    """

    script:
    """
    cp ${params.exon_bed} exons.bed
    ${HAL}/bin/hal4dExtract --hdf5InMemory ${params.hal} ${params.reference} exons.bed neutralRegions.bed 
    """
}

process phyloFit {
    tag "pFit"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path input
    path HAL
    //path PHAST

    output:
    path "neutralModel.mod"
    
    
    stub:
    """
    touch neutralModel.mod
    """

    script:
    if (params.hal4d)
    """
    mytree=`${HAL}/bin/halStats --tree ${params.hal} | python -c "import sys, re; pattern = r'Inner[0-9]*'; sys.stdout.write( re.sub(r':[0-9].[0-9]*', '', re.sub(pattern,'', sys.stdin.readline()).replace(':0)Anc00',''))[1:] ) "`
    phyloFit --tree \$mytree --msa-format SS --subst-mod SSREV --sym-freqs --precision HIGH --out-root neutralModel ${input} 
    """
    else
    """
    tree=`${HAL}/bin/halStats --tree ${params.hal}`
    phyloFit --out-root neutralModel --tree "\$tree" ${input}
    """
}

process makeMaf {
    tag "maf"
    publishDir "${params.outdir}/PHYLOP/MAF/CHROM", mode: "${params.publish_dir_mode}", overwrite: true
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path HAL
    tuple val(N), val(chrom)

    output:
    path "alignments_${chrom}.maf"
    
    
    stub:
    """
    touch alignments_${chrom}.maf
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    ${HAL}/bin/hal2maf --refSequence ${chrom} \
        --refGenome ${params.reference} --noAncestors --noDupes \
        --hdf5InMemory ${params.hal} alignments_${chrom}.maf
    """
}

process combine_mafs {
    tag "combmaf"
    publishDir "${params.outdir}/PHYLOP/MAF", mode: "${params.publish_dir_mode}", overwrite: true
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path mafs

    output:
    path "alignments.maf"
    
    
    stub:
    """
    touch alignments.maf
    """

    script:
    $/
    import sys, os
    infiles = [i for i in os.path.listdir() if '.maf' in i]
    outfile = open('alignments.maf', 'w')
    for n, infile in enumerate(infiles):
        for line in open(infile):
            if '#' in line[0] and n != 0:
                continue
            outfile.write(line)
    /$
}


process make4dmaf {
    tag "maf"
    publishDir "${params.outdir}/PHYLOP/MAF", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_vlargemem"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path bedfile
    path HAL

    output:
    path "4d.maf"
    
    stub:
    """
    touch 4d.maf
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    ${HAL}/bin/hal2mafMP.py \
        --noDupes \
        --noAncestors \
        --numProc ${task.cpus} \
        --refGenome ${params.reference} \
        --refTargets ${bedfile} \
        --hdf5InMemory ${params.hal} 4d.maf
    sed -i -e 2d 4d.maf
    """
}

process msa_view {
    tag "msa"
    publishDir "${params.outdir}/PHYLOP/MSA", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_largemem"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path maf
    path HAL

    output:
    path "${maf.simpleName}.ss"
    
    
    stub:
    """
    touch ${maf.simpleName}.ss
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    genomes=`${HAL}/bin/halStats ${params.hal} --genomes | python -c 'import sys; a = [i for i in sys.stdin.readline().split() if "Inner" not in i and "Anc" not in i]; print(",".join(a))'`
    msa_view -o SS -z --in-format MAF --aggregate \$genomes ${maf} > ${maf.simpleName}.ss
    """
}


process phastCons {
    tag "pcon"
    publishDir "${params.outdir}/PHYLOP/PHASTCONS/CHR", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_multi"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path HAL
    path maf
    path model

    output:
    path "${maf.simpleName}.cons.bed"
    path "${maf.simpleName}.scores.wig"
    
    
    stub:
    """
    touch ${maf.simpleName}.cons.bed
    touch ${maf.simpleName}.scores.wig
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    phastCons --most-conserved ${maf.simpleName}.cons.bed ${maf} ${model} > ${maf.simpleName}.scores.wig
    """
}

process collect{
    tag "col"
    publishDir "${params.outdir}/PHYLOP/PHASTCONS", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_multi"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path beds

    output:
    path "conserved.bed"
    
    
    stub:
    """
    touch conserved.bed
    """

    script:
    """
    cat ${beds} | sed 's/alignments_//g' | bedtools sort -i - > conserved.bed
    """
}


process phyloPtrain {
    tag "pptrain"
    publishDir "${params.outdir}/PHYLOP/NEUTRAL", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium_multi"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path HAL
    path ss

    output:
    path "neutralModel.mod"
    
    
    stub:
    """
    touch neutralModel.mod
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    ${HAL}/bin/halPhyloPTrain.py ${params.hal} ${params.reference} ${neutral} neutralModel.mod --numProc ${task.cpus} 
    """
}


process phyloP {
    tag "ppmp"
    publishDir "${params.outdir}/PHYLOP/BED/CHR", mode: "${params.publish_dir_mode}", overwrite: true
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path HAL
    path model
    tuple val(n), val(chr)

    output:
    path "constrained_elements_${chr}.bed"
    
    
    stub:
    """
    touch constrained_elements_${chr}.bed
    """

    script:
    """
    export PYTHONPATH=${HAL}/lib:${HAL}/submodules/sonLib/src
    export PATH=\$PATH:${HAL}/bin
    ${HAL}/bin/halPhyloPMP.py ${params.hal} ${params.reference} ${model} constrained_elements_${chr}.wig \
        --refSequence ${chr} \
        --chromSizes ${params.reference}.sizes \
        --numProc ${task.cpus} 
    wigToBigWig constrained_elements_${chr}.wig ${params.reference}.sizes /dev/stdout | \
        bigWigToBedGraph /dev/stdin /dev/stdout | \
        bedtools sort -i - > constrained_elements_${chr}.bed && \
        rm constrained_elements_${chr}.wig
    """
}


process bigWigToBed {
    tag "bw2bed"
    publishDir "${params.outdir}/PHYLOP/BED/CHR", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda "${baseDir}/envs/phast_environment.yml"

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
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path beds

    output:
    path "constrained.bed"
    
    
    stub:
    """
    touch constrained.bed
    """

    script:
    """
    cat ${beds} | bedtools sort -i - > constrained.bed
    """
}


process get_constrained {
    tag "constr"
    publishDir "${params.outdir}/PHYLOP/CONSTRAINED", mode: "${params.publish_dir_mode}", overwrite: true
    label "largemem"
    conda "${baseDir}/envs/phast_environment.yml"

    input:
    path bed

    output:
    path "constrained.bed"
    
    
    stub:
    """
    touch constrained.bed
    """

    script:
    """
    cat ${bed} | bedtools sort -i - > constrained.bed
    """
}


process filter {
    tag "filt"
    publishDir "${params.outdir}/PHYLOP/VCF", mode: "${params.publish_dir_mode}", overwrite: true
    label "medium"

    input:
    path vcf
    path tbi
    path bed

    output:
    path "${vcf.simpleName}.unconstrained.vcf.gz"
    path "${vcf.simpleName}.unconstrained.vcf.gz.tbi"
    
    
    stub:
    """
    touch ${vcf.simpleName}.unconstrained.vcf.gz
    touch ${vcf.simpleName}.unconstrained.vcf.gz.tbi
    """

    script:
    """
    bedtools intersect -header -v -a ${vcf} -b ${bed} | bgzip -c > ${vcf.simpleName}.unconstrained.vcf.gz
    tabix -p vcf ${vcf.simpleName}.unconstrained.vcf.gz
    """
}

