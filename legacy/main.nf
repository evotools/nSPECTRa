#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 * if Nsnps is set to 0, use all
 */
params.variants = "file.vcf.gz"
params.hal = 'bovinae.hal'
params.reference = 'hereford'
params.target = "Anc3"
params.outdir = 'COW'
params.sif = 'PATH/TO/hal.sif'
params.annotation = ''


/*
 * Phase 1. From hal to ancestral fasta  
 */
process hal2maf {
    tag "hal2maf"
    publishDir "./${params.outdir}/MAF", mode: 'symlink', overwrite: true

    memory { 64.GB * task.attempt }
    time { 12.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"
    
    output:
    path "${params.reference}_${params.target}.maf"  
    
    script:
    """
    ln ${params.hal} 
    halname=`basename ${params.hal}`
    singularity exec --bind $PWD:/mnt ${params.sif} hal2maf \
        --refGenome ${params.reference} \
        --targetGenomes ${params.target} \
        --hdf5InMemory \$halname ${params.reference}_${params.target}.maf
    """
}

process maf2fasta {
    tag "maf2fasta"
    publishDir "./${params.outdir}/ANCFASTA", mode: 'symlink', overwrite: true

    memory { 16.GB * task.attempt }
    time { 6.hour * task.attempt }
    cpus { 4 * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"
    module "roslin/samtools/1.10"

    beforeScript = """
    . /etc/profile.d/modules.sh
    sleep 2;
    module load anaconda
    source activate mutyperenv
    """

    input:
    path maf

    output:
    path "${params.reference}_${params.target}.fasta", emit: ancfasta
    path "${params.reference}_${params.target}.fasta.fai", emit: ancfai
    path "${params.reference}_${params.target}_regions.bed", emit: ancbed
    
    script:
    """
    MAF2ANCFA --maf ${maf} \
        --threads ${task.cpus} \
        --out ${params.reference}_${params.target}
    samtools faidx ${params.reference}_${params.target}.fasta
    awk 'BEGIN{OFS="\t"};{print \$1,"0",\$2}' > ${params.reference}_${params.target}_regions.bed
    """
}

/*
 * Phase 2: Change and lift positions
 */
process liftToAncestor {
    tag "liftover"

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    output:
    path "${params.reference}_variantslifted_1nt.bed", emit: varlifted
    path "${params.reference}_MISSED_1nt.txt", emit: varmissed

    script:
    """
    ln ${params.hal} 
    bcftools view -v snps -m2 -M2 ${params.variants} | sed 's/_pilon//g' |\
        awk -v val=0 'BEGIN{OFS="\t"}; \$1!~"#" {print \$1,\$2-1-val,\$2+val,\$1"_"\$2"_"\$4"_"\$5}' | \
        sed 's/,/-/g' > ${params.reference}_variants2lift_1nt.bed
    halname=`basename ${params.hal}`
    singularity exec --bind $PWD:/tmp ${params.sif} halLiftover \
        --hdf5InMemory \$halname \
        ${params.reference} ${params.reference}_variants2lift_1nt.bed \
        ${params.target} ${params.reference}_variantslifted_1nt.bed 2> ${params.reference}_MISSED_1nt.txt
    """
}


/*
 * Phase 3: Generate ancestral genome to process
 */
process makeRefTgtFasta {
    tag "reftgtfasta"

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    output:
    path "${params.reference}.fasta", emit: reffasta
    path "${params.target}.fasta", emit: tgtfasta
    path "${params.reference}.fasta", emit: reffai
    path "${params.target}.fasta", emit: tgtfai

    script:
    """
    ln ${params.hal} 
    halname=`basename ${params.hal}`
    singularity exec --bind $PWD:/mnt ${params.sif} hal2fasta \
        --hdf5InMemory \$halname ${params.target} > ${params.target}.fasta
    singularity exec --bind $PWD:/mnt ${params.sif} hal2fasta \
        --hdf5InMemory \$halname ${params.reference} > ${params.reference}.fasta
    samtools faidx ${params.target}.fasta
    samtools faidx ${params.reference}.fasta
    """
}

/*
 * Phase 5: make annotation for variants 
 */
process makeAnnotation {
    tag "makeAnnot"

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    input:
    path liftedbed
    path tgtfasta

    output:
    path "Ancestral_annotation_${params.reference}.txt.gz", emit: annotation
    path "Ancestral_annotation_${params.reference}.txt.gz.tbi", emit: annottbi

    script:
    """
    bedtools getfasta -fi ${tgtfasta} -bed ${liftedbed} -name -tab |
        GetRightAllele stdin | \
        awk 'BEGIN{OFS="\t"}; NR>1{print \$1,\$2,\$4}' | \
        bgzip -c > Ancestral_annotation_${params.reference}.txt.gz && \
        tabix -b 2 -e 2 -s 1 Ancestral_annotation_${params.reference}.txt.gz
    """

}

/*
 * Phase 6: variants filtering
 */
process filtering {
    tag "filter"

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    cpus { 4 * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    output:
    path "missingness_het_${params.reference}.txt"

    script:
    """
    if [[ ${params.variants} == *.bcf ]]; then
        plink --cow --bcf ${params.variants} --allow-extra-chr --mind 0.1 --double-id --het --out missingness_het_${params.reference} --threads ${task.cpus}
    else
        plink --cow --vcf ${params.variants} --allow-extra-chr --mind 0.1 --double-id --het --out missingness_het_${params.reference} --threads ${task.cpus}
    fi
    awk 'NR>1 && ((\$5-\$3) / \$5) > 0.15 && ((\$5-\$3) / \$5 )<0.85 {print \$2}' missingness_het_${params.reference}.het > missingness_het_${params.reference}.txt
    """
}

/*
 * Phase 7: run mutyper
 */
process mutyper {
    tag "mutyper"
    publishDir "chromosomeSpectra", mode: 'copy', overwrite: true

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    input:
    path filter
    path annotation
    path annotationtbi
    path ancfasta
    path ancfai
    val k

    output:
    path "mutationSpectra_${params.reference}*.txt"

    script:
    if (params.annotation)
    """
    echo "Run mutyper"
    source activate mutyperenv
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > header_${params.reference}.txt
    bcftools view -v snps -m2 -M2 -S ${filter} ${params.variants} |
        sed 's/_pilon//g' |
        bcftools annotate -a ${annotation} -h header_${params.reference}.txt -c CHROM,POS,AA |
        vcftools --vcf - --minGQ 30 --max-missing 0.90 --recode --recode-INFO-all --stdout |
        mutyper variants -k ${k} --strand_file ${params.annot} ${ancfasta} - |
        mutyper spectra - > mutationSpectra_${params.reference}_${region}.txt
    """
    else
    """
    echo "Run mutyper"
    source activate mutyperenv
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > header_${params.reference}.txt
    bcftools view -v snps -m2 -M2 -S ${filter} ${params.variants} |
        sed 's/_pilon//g' |
        bcftools annotate -a ${annotation} -h header_${params.reference}.txt -c CHROM,POS,AA |
        vcftools --vcf - --minGQ 30 --max-missing 0.90 --recode --recode-INFO-all --stdout |
        mutyper variants -k ${k} ${ancfasta} - |
        mutyper spectra - > mutationSpectra_${params.reference}_${region}.txt
    """
}

process group_results {
    tag "group_results"
    publishDir "mutationSpectra", mode: 'copy', overwrite: true

    memory { 8.GB * task.attempt }
    time { 6.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"
    modules "R/3.5.3"

    input:
    file data

    output:
    file "mutationSpectra_${params.reference}.txt"

    script:
    $/
    #!/usr/bin/env Rscript
    allmuts = list.files(pattern = ".txt", recursive = T)
    if (exists("mutSpectra")){rm(mutSpectra)}
    for (f in allmuts){
        mytab = read.table(f, h=T)
        if (!exists("mutSpectra")){
            mutSpectra = mytab
        } else {
            mutSpectra = merge(mutSpectra, mytab, all.x = T, all.y = T)
        }
    }
    write.table(mutSpectra, "mutationSpectra_${params.reference}.txt", sep = ",", col.names = T, row.names = F, quotes = F)
    /$
}


/*
 * Start workflow
 */
workflow {
    /*Create ancestral fasta file*/
    hal2maf()
    maf2fasta(hal2maf.out)

    /* Create lifted positions */
    liftToAncestor()
    makeRefTgtFasta()

    /* Make annotation file for the variants*/
    makeAnnotation( liftToAncestor.out.varlifted, makeRefTgtFasta.out.tgtfasta )
    filtering()

    /* Start ,utyper on each chromosome separately */
    kvals = Channel
        .from( 3, 5, 7 )
    mutyper( filtering.out, makeAnnotation.out.annotation, makeAnnotation.out.annottbi, maf2fasta.out.ancfasta, maf2fasta.out.ancfai, kvals )

    /* Collect outputs */
    group_results(mutyper.out.collect())

}