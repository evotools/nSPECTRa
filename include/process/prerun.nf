
// Get chromosome list
process chromosomeList {
    tag "chr"
    label "medium"

    input:
    path ch_vcf
    path ch_tbi

    output:
    path "sequences.csv"

    
    stub:
    """
    echo "1,seq1" > sequences.csv
    echo "2,seq2" >> sequences.csv
    echo "3,seq3" >> sequences.csv
    """

    script:
    """
    tabix -l ${ch_vcf} | awk 'BEGIN{OFS=","};{print NR, \$0}' > sequences.csv 
    """
}

// Exzclude lists of small populations
process get_masks{
    label "small"

    input:
    path genome

    output:
    path "masks.bed"

    stub:
    """
    touch masks.bed
    """

    script:
    """
    faToTwoBit ${genome} ${genome.simpleName}.2bit
    twoBitInfo -maskBed ${genome.simpleName}.2bit masks.bed
    """

}


// Split vcf by chromosome
process splitvcf {
    tag "split"
    label "small"

    input:
    tuple path(ch_vcf), path(ch_tbi)
    tuple val(index), val(chrom)

    output:
    tuple val(chrom), path("subset_${chrom}.vcf.gz")

    
    stub:
    """
    touch subset_${chrom}.vcf.gz
    touch subset_${chrom}.vcf.gz.tbi
    """

    script:
    """
    bcftools view --threads ${task.cpus} -r ${chrom} -O z -o subset_${chrom}.vcf.gz ${ch_vcf} && \
        tabix -p vcf subset_${chrom}.vcf.gz
    """
}

process get_individuals {
    tag "inds"
    label "small"

    input:
    path ch_vcf
    path ch_tbi

    output:
    path "individuals.txt"
    
    
    stub:
    """
    echo "FID1_IID1" > individuals.txt
    echo "FID1_IID2" >> individuals.txt
    echo "FID2_IID1" >> individuals.txt
    echo "FID2_IID2" >> individuals.txt
    """

    script:
    """
    vcfsamplenames ${ch_vcf} > individuals.txt
    """
}


// Get breeds' names
process get_breeds {
    tag "breeds"
    label "small"

    input:
    path indvs

    output:
    path "breeds.csv"
    
    
    stub:
    """
    awk 'BEGIN{OFS=","; FS="_"}; {print \$1}' ${indvs} | sort | uniq | awk '{print \$NR, \$0}; END{print \$NR+1, "ALL"}' > breeds.csv 
    """

    script:
    """
    awk 'BEGIN{OFS=","; FS="_"}; {print \$1}' ${indvs} | sort | uniq | awk '{print \$NR, \$0}; END{print \$NR+1, "ALL"}' > breeds.csv 
    """
}

// Run shapeit2/beagle for each chromosome
process beagle {
    tag "beagle.${chrom}"
    label "large_largemem"

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path beagle

    output:
    tuple val(chrom), path("prephase_${chrom}.vcf.gz"), path("prephase_${chrom}.vcf.gz.tbi")

    
    stub:
    """
    touch prephase_${chrom}.vcf.gz
    touch prephase_${chrom}.vcf.gz.tbi
    """

    script:
    def ne = params.neval ? "ne=${params.neval}" : ""
    """
    javamem=`python -c "import sys; maxmem=int(sys.argv[1]); print( maxmem - int(maxmem * .1) )" ${task.memory.toGiga()}`
    java -jar -Xmx\${javamem}G ${beagle} gt=${vcf} ${ne} chrom=${chrom} out=prephase_${chrom} nthreads=${task.cpus}
    if [ ! -e prephase_${chrom}.vcf.gz.tbi ]; then
        tabix -p vcf prephase_${chrom}.vcf.gz
    fi
    """
}

process shapeit4 {
    tag "impute.${chrom}"
    label "medium_multi"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("prephase_${chrom}.vcf.gz"), path("prephase_${chrom}.vcf.gz.tbi")

    
    stub:
    """
    touch prephase_${chrom}.vcf.gz
    touch prephase_${chrom}.vcf.gz.tbi
    """

    script:
    def psfield = params.whatshap ? "--use-PS 0.0001" : null
    def shapeit = params.shapeit ? "${params.shapeit}" : "shapeit"
    def ne = params.neval ? "--effective-size ${params.neval}" : ""
    """
    ${shapeit} --input ${vcf} -T ${task.cpus} --region ${chrom} ${ne} --output prephase_${chrom}.vcf.gz --sequencing ${psfield}
    if [ ! -e prephase_${chrom}.vcf.gz.tbi ]; then
        tabix -p vcf prephase_${chrom}.vcf.gz
    fi
    """
    // if (params.shapeit)
    // """
    // ${params.shapeit} --input ${vcf} -T ${task.cpus} --region ${chrom} --effective-size ${params.neval} --output prephase_${chrom}.vcf.gz --sequencing ${psfield}
    // """
    // else
    // """
    // shapeit4 --input ${vcf} -T ${task.cpus} --region ${chrom} --effective-size ${params.neval} --output prephase_${chrom}.vcf.gz --sequencing ${psfield}
    // """
}

process shapeit5 {
    tag "impute.${chrom}"
    label "medium_multi"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("prephase_${chrom}.vcf.gz"), path("prephase_${chrom}.vcf.gz.tbi")

    
    stub:
    """
    touch prephase_${chrom}.vcf.gz
    touch prephase_${chrom}.vcf.gz.tbi
    """

    script:
    def shapeit = params.shapeit ? "${params.shapeit}" : "phase_common"
    def ne = params.neval ? "--hmm-ne ${params.neval}" : ""
    """
    ${shapeit} --input ${vcf} --region ${chrom} --filter-maf 0.001 --output prephase_${chrom}.vcf.gz --thread ${task.cpus}
    """
}

// Split VCf by chromosome
process split_vcf {
    tag "split.${chrom}"
    label "small"

    input:
    tuple val(idx), val(chrom)
    path vcf
    path tbi

    output:
    tuple val(chrom), path("prephase_${chrom}.vcf.gz")

    
    stub:
    """
    touch prephase_${chrom}.vcf.gz
    """

    script:
    """
    bcftools view --threads ${task.cpus} -O z -r ${chrom} ${vcf} > prephase_${chrom}.vcf.gz
    """
}

process combineVcf {
    tag "combine"
    label 'medium'
    publishDir "${params.outdir}/imputed", mode: "${params.publish_dir_mode}", overwrite: true


    input:
        path "vcfs/*"

    output:
        path "IMPUTED.vcf.gz" 
        path "IMPUTED.vcf.gz.tbi"
  
    
    stub:
    """
    touch IMPUTED.vcf.gz
    touch IMPUTED.vcf.gz.tbi
    """

    script:
    """
    bcftools concat -O u vcfs/*.vcf.gz | bcftools sort -T ./ -O z -m 5G > IMPUTED.vcf.gz && \
        tabix -p vcf IMPUTED.vcf.gz
    """
}


/*
 * Step 5. Get identity by descent through beagle
 */
process ibd {
    tag "ibd"
    publishDir "${params.outdir}/ibd", mode: "${params.publish_dir_mode}", overwrite: true


    cpus 4
    memory { 16.GB * task.attempt }
    time { 24.hour * task.attempt }
    clusterOptions "-P roslin_ctlgh -l h_vmem=${task.memory.toString().replaceAll(/[\sB]/,'')}"

    input: 
        tuple val(chrom), path(imputed), path(imputed_tbi)

    output: 
        tuple val(chrom), path("IBD.${chrom}.ibd.gz")
        tuple val(chrom), path("IBD.${chrom}.hbd.gz")
  
    
    stub:
    """
    touch IBD.${chrom}.ibd.gz
    touch IBD.${chrom}.hbd.gz
    """

    script:
    """
    # See 10.1534/genetics.113.150029 for details on filtering
    javamem=`python -c "import sys; maxmem=int(sys.argv[1]); print( maxmem - int(maxmem * .1) )" ${task.memory.toGiga()}`
    java -Xmx\${javamem}G -jar ${params.refinedibd} \
            gt=${imputed} \
            out=IBD.${chrom} \
            nthreads=${task.cpus} \
            ${params.ibd_params}
    """
}

process collect_vcf {

    input:
    path path(subset)

    output:
    path "ready.vcf.gz"

    
    stub:
    """
    touch ready.vcf.gz
    """

    script:
    """
    vcf-concat ${vep_out} | vcf-sort | bgzip -c > ready.vcf.gz 
    """

}

// Make shapeit files
process make_shapeit{
    tag "shapeit.${chrom}"

    input:
    tuple val(chrom), path(vcf)

    output:
    tuple val(chrom), path("${chrom}.SHAPEIT.*")

    
    stub:
    """
    touch ${chrom}.SHAPEIT.T
    """

    script:
    """
    vcftools --gzvcf ${vcf} --IMPUTE --out ${chrom}.SHAPEIT
    """
}


// Compute derived allele freq.
process daf {
    tag "mutyper"
    label "medium"
    publishDir "${params.outdir}/mutyper/daf", mode: "${params.publish_dir_mode}", overwrite: true


    input:
    path vcf
    path tbi

    output:
    path "daf.csv.gz"

    
    stub:
    """
    echo "" | bgzip -c > daf.csv.gz
    """

    script:
    """
    bcftools +fill-tags ${vcf} -- -t AF,AN,AC | bcftools query -H -f "%CHROM,%POS,%AF\n" | bgzip -c > daf.csv.gz
    """
}

// Smile plot for the derived allele frequencies
process smile {
    label "renv"
    publishDir "${params.outdir}/mutyper/smile", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path daf

    output:
    path "smile.pdf"
    path "smile.png"

    stub:
    """
    touch smile.pdf
    touch smile.png
    """

    script:
    """
    #!/usr/bin/env Rscript
    if (TRUE){
        library(tidyverse)
        library(ggpubr)
        library(ggsci)
    }
    setwd('~/Analysis/DAF/')
    daf <- read_csv('${daf}', col_names = c('chrom', 'pos', 'af'), comment = '#')

    p <- gghistogram(daf, x="af", bins=100, color = "#00AFBB", fill = "#00AFBB")
    ggsave('smile.pdf', device = 'pdf', heigth=9, width=16)
    ggsave('smile.png', device = 'png', heigth=9, width=16)
    """
}