
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
    cpus = 4
    memory = { params.greedy ? 96.GB * task.attempt : 32.GB * task.attempt }
    time = { 23.h * task.attempt }

    input:
    path bed
    path fasta
    path fai
    tuple val(contig), val(header)

    output:
    tuple path("./${contig}_ancestral_states.bed.gz"), path("${contig}.fasta"), optional: true
    
    
    stub:
    """
    touch ./${contig}_ancestral_states.bed.gz
    touch ${contig}.fasta
    """

    script:
    def greedy = params.greedy ? "--greedy" : ""
    """
    samtools faidx ${fasta} ${contig} > ${contig}.fasta
    VERTICALIZE -b ${bed} -t ${task.cpus} --region ${contig} -f ${contig}.fasta -o /dev/stdout |\
        CONSENSE -b /dev/stdin -t ${task.cpus} --region ${contig} -f ${contig}.fasta -o ${contig}_ancestral_states.bed ${greedy}
    bgzip ${contig}_ancestral_states.bed
    """
}

/*
 * Phase 3: Generate ancestral genome to process
 */
process makeRefTgtFasta {
    tag "reftgtfasta"
    label "medium"
    container { params.cactus_version ? "quay.io/comparative-genomics-toolkit/cactus:${params.cactus_version}" : "quay.io/comparative-genomics-toolkit/cactus:latest" }

    input:
    path HAL

    output:
    path "${params.reference}.fasta", emit: reffasta
    path "${params.target}.fasta", emit: tgtfasta
    path "${params.reference}.fasta.fai", emit: reffai
    path "${params.target}.fasta.fai", emit: tgtfai

    
    stub:
    """
    echo '>seq1.r' > ${params.reference}.fasta
    echo 'AACCTTGG' >> ${params.reference}.fasta
    echo '>seq2.r' >> ${params.reference}.fasta
    echo 'AACCTTGT' >> ${params.reference}.fasta
    echo '>seq3.r' >> ${params.reference}.fasta
    echo 'TTCCTTGT' >> ${params.reference}.fasta
    samtools faidx ${params.reference}.fasta

    echo '>seq1.t' > ${params.target}.fasta
    echo 'AACCTTGG' >> ${params.target}.fasta
    echo '>seq2.t' >> ${params.target}.fasta
    echo 'AACCTTGT' >> ${params.target}.fasta
    echo '>seq3.t' >> ${params.target}.fasta
    echo 'TTCCTTGT' >> ${params.target}.fasta
    samtools faidx ${params.target}.fasta
    """

    script:
    """
    hal2fasta --hdf5InMemory ${HAL} ${params.target} > ${params.target}.fasta
    hal2fasta --hdf5InMemory ${HAL} ${params.reference} > ${params.reference}.fasta
    samtools faidx ${params.target}.fasta
    samtools faidx ${params.reference}.fasta
    """
}

process filter_by_size {
    tag "reffilt"
    label "medium"

    input:
    path REF
    path FAI

    output:
    path "${params.reference}.large.fasta", emit: reffasta
    path "${params.reference}.large.fasta.fai", emit: reffai
    
    
    stub:
    """
    touch ${params.reference}.large.fasta
    touch ${params.reference}.large.fasta.fai
    """

    script:
    """
    samtools faidx ${REF} \$( awk -v size=${params.ref_min_size} '\$2>=size {print \$1}' ${FAI} ) > ${params.reference}.large.fasta && 
        samtools faidx ${params.reference}.large.fasta
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
    path "${reffa.baseName}.anc.fa"

    
    stub:
    """
    touch ${reffa.baseName}.anc.fa
    """

    script:
    """
    BED2ANCFASTA -b ${bedfile} -f ${reffa} -o ${reffa.baseName}.anc.fa
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

// Rename fasta file if requested
process rename_fasta {
    label "small"
    publishDir "${params.outdir}/ancestral/ancestral_genome_renamed", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path fasta
    path conversion_table

    output:
    path "${fasta.simpleName}.renamed.fasta"
    path "${fasta.simpleName}.renamed.fasta.fai"

    
    stub:
    """
    touch ${fasta.simpleName}.renamed.fasta
    touch ${fasta.simpleName}.renamed.fasta.fai
    """

    script:
    """
    FASTA_REHEADER ${fasta} ${conversion_table} > ${fasta.simpleName}.renamed.fasta
    samtools faidx ${fasta.simpleName}.renamed.fasta
    """
}