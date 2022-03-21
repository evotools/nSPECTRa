
// process maf2fasta {
//     tag "maf2fasta"
//     publishDir "${params.outdir}/ancestral/ancfasta", mode: 'symlink', overwrite: true

//     memory { 16.GB * task.attempt }
//     time { 6.hour * task.attempt }
//     cpus { 4 * task.attempt }
//     clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"
//     module "roslin/samtools/1.10"

//     beforeScript = """
//     . /etc/profile.d/modules.sh
//     sleep 2;
//     module load anaconda
//     source activate mutyperenv
//     """

//     input:
//     path maf

//     output:
//     path "${params.reference}_${params.target}.fasta", emit: ancfasta
//     path "${params.reference}_${params.target}.fasta.fai", emit: ancfai
//     path "${params.reference}_${params.target}_regions.bed", emit: ancbed
    
//     script:
//     """
//     MAF2ANCFA --maf ${maf} \
//         --threads ${task.cpus} \
//         --out ${params.reference}_${params.target}
//     samtools faidx ${params.reference}_${params.target}.fasta
//     awk 'BEGIN{OFS="\t"};{print \$1,"0",\$2}' > ${params.reference}_${params.target}_regions.bed
//     """
// }

process maf2bed {
    label "large"
    publishDir "${params.outdir}/ancestral/ancbed", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    path maf

    output:
    path "${params.reference}_${params.target}.bed"
    
    
    stub:
    """
    touch ${params.reference}_${params.target}.bed
    """

    script:
    """
    MAFSPLIT -m ${maf} -o stdout > tmp.split.maf
    mafSTRANDER --maf tmp.split.maf --strand + > ${params.reference}_${params.target}.split.plus.maf && rm tmp.split.maf
    MAF2ANCFA -m ${params.reference}_${params.target}.split.plus.maf -d ${params.reference} -t ${task.cpus} -O bed -o ${params.reference}_${params.target}
    """
}

process bed2vbed{
    label "large"

    input:
    path bed
    path contig

    output:
    tuple path("./${contig.simpleName}_ancestral_states.bed.gz"), path(contig)
    
    
    stub:
    """
    touch ./${contig.simpleName}_ancestral_states.bed.gz
    """

    script:
    """
    awk -v chrid=${contig.simpleName} '\$1==chrid {print}' ${bed} > ./${contig.simpleName}.bed
    BED2VBED -b ./${contig.simpleName}.bed | \
            sort --buffer-size=25G --parallel=${task.cpus} -T ./TMP/ -k1,1 -k2,2n - | \
            COMBINE -b - -f ${contig} | \
            CONSENSE -b - | \
            bgzip -c > ${contig.simpleName}_ancestral_states.bed.gz && rm ./${contig.simpleName}.bed
    """
}

/*
 * Phase 3: Generate ancestral genome to process
 */
process makeRefTgtFasta {
    tag "reftgtfasta"
    label "medium"

    input:
    path HAL
    path CACTUS

    output:
    path "${params.reference}.fasta", emit: reffasta
    path "${params.target}.fasta", emit: tgtfasta
    path "${params.reference}.fasta.fai", emit: reffai
    path "${params.target}.fasta.fai", emit: tgtfai

    
    stub:
    """
    touch ${params.reference}.fasta
    touch ${params.target}.fasta
    touch ${params.reference}.fasta.fai
    touch ${params.target}.fasta.fai
    """

    script:
    """
    ${CACTUS}/bin/hal2fasta --hdf5InMemory ${HAL} ${params.target} > ${params.target}.fasta
    ${CACTUS}/bin/hal2fasta --hdf5InMemory ${HAL} ${params.reference} > ${params.reference}.fasta
    samtools faidx ${params.target}.fasta
    samtools faidx ${params.reference}.fasta
    """
}

process splitfasta {
    label "medium"

    input:
    path fasta
    path fai

    output:
    path "*.sub.fa"

    
    stub:
    """
    touch seq1.sub.fa
    touch seq2.sub.fa
    touch seq3.sub.fa
    """

    script:
    """
    for i in `awk '{print \$1}' ${fai}`; do
        samtools faidx ${fasta} \$i > ./\$i.sub.fa
    done
    """
}


process bed2ancfa {
    label "medium"

    input:
    tuple path(bedfile), path(reffa)

    output:
    path "${reffa.simpleName}.anc.fa"

    
    stub:
    """
    touch ${reffa.simpleName}.anc.fa
    """

    script:
    """
    BED2ANCFASTA -b ${bedfile} -f ${reffa} -o ${reffa.simpleName}.anc.fa
    """
}

process collectAncfa {
    label "small"
    publishDir "${params.outdir}/ancestral/ancestral_genome", mode: "${params.publish_dir_mode}", overwrite: true



    input:
    path ancfas

    output:
    path "ancestral.fasta"
    path "ancestral.fasta.fai"

    
    stub:
    """
    touch ancestral.fasta
    touch ancestral.fasta.fai
    """

    script:
    """
    cat ${ancfas} > ancestral.fasta
    samtools faidx ancestral.fasta
    """
}

process makefai {
    label "small"
    publishDir "${params.outdir}/ancestral/ancestral_genome", mode: "${params.publish_dir_mode}", overwrite: true



    input:
    path ancfa

    output:
    path "*.fai"

    
    stub:
    """
    touch ${ancfa.simpleName}.fai
    """

    script:
    """
    samtools faidx ${ancfa}
    """
}