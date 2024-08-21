/*
 * Phase 2: Change and lift positions
 */
process liftToAncestor {
    tag "liftover"
    label "medium"

    output:
    path "${params.reference}_variantslifted_1nt.bed", emit: varlifted
    path "${params.reference}_MISSED_1nt.txt", emit: varmissed

    
    stub:
    """
    touch ${params.reference}_variantslifted_1nt.bed
    touch ${params.reference}_MISSED_1nt.txt
    """

    script:
    """
    ln ${params.hal} 
    bcftools view --threads ${task.cpus} -v snps -m2 -M2 ${params.variants} | sed 's/_pilon//g' |\
        awk -v val=0 'BEGIN{OFS="\t"}; \$1!~"#" {print \$1,\$2-1-val,\$2+val,\$1"_"\$2"_"\$4"_"\$5}' | \
        sed 's/,/-/g' > ${params.reference}_variants2lift_1nt.bed
    halname=`basename ${params.hal}`
    singularity exec --bind $PWD:/tmp ${params.sif} halLiftover \
        --hdf5InMemory \$halname \
        ${params.reference} ${params.reference}_variants2lift_1nt.bed \
        ${params.target} ${params.reference}_variantslifted_1nt.bed 2> ${params.reference}_MISSED_1nt.txt
    """
}



process makeAncestralSequence {
    tag "makeAncSeq"
    label "medium"

    input:
    path vcf
    path reference
    path target
    path chain

    output:
    path "ancestral_sequence.fasta"
    path "ancestral_sequence.fasta.fai"

    
    stub:
    """
    touch ancestral_sequence.fasta
    touch ancestral_sequence.fasta.fai
    """

    script:
    if (params.converter == "mutyper")
    """
    mutyper ancestor ${vcf} ${reference} ${outgroup} ${chain} ancestral_sequence.fasta
    samtools faidx ancestral_sequence.fasta
    """
    else
    """
    module load java
    git clone https://github.com/abyzovlab/vcf2diploid
    cd vcf2diploid && make && ..
    java -jar ./vcf2diploid/vcf2diploid.jar -id ${params.target} -chr ${reference} -vcf ${vcf}
    """
}
