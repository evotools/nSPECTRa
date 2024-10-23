
/* 
 * Calculate the Fst. 
 * Use bcftools speed to extract the chromosome of interest, 
 * pipe it to vcftools and calculate the FST
 */
// Prepare gone inputs
process gone_inputs { 
    publishDir "${params.outdir}/ne/input/${pop}", mode: 'copy', overwrite: true
    label "medium"

    input: 
        path vcf
        path tbi
        tuple val(pop), val(samplefile)

    output: 
        tuple val(pop), path("${pop}.ped"), path("${pop}.map")
  
    stub:
    """
    touch ${pop}.ped
    touch ${pop}.map    
    """

    script:
    """
    awk -v pop=${pop} '{print pop, \$1}' ${samplefile} > keep.txt
    plink --allow-extra-chr --chr-set 90 --vcf ${vcf} --keep keep.txt --const-fid ${pop} --recode --out ${pop} --maf 0.05 --geno 0.05 --threads ${task.cpus}
    mv ${pop}.map ${pop}.bck
    FixMap ${pop}.bck ${pop}
    """
}


process gone_run { 
    publishDir "${params.outdir}/ne/results/${pop}", mode: 'copy', overwrite: true
    label "large"

    input: 
        tuple val(pop), path(ped), path(map)
        path GONE

    output: 
        path "${pop}.ne"
  
    stub:
    """
    touch ${pop}.ne
    """

    script:
    """
    min_by_chr=`cut -f 1 ${pop}.map | sort | uniq -c | sort -k1,1n`

    # Get GONE
    cat ${GONE}/Linux/INPUT_PARAMETERS_FILE | \
        sed "s/threads=-99/threads=${task.cpus}/g" | \
        sed "s/maxNSNP=-99/maxNSNP=${params.ne_subset}/g" > INPUT_PARAMETERS_FILE

    if [ ! -e script_GONE.sh ]; then
        cp ${GONE}/Linux/script_GONE.sh ./
        chmod a+x ./script_GONE.sh
    fi

    if [ ! -e PROGRAMMES ]; then
        cp -r ${GONE}/Linux/PROGRAMMES ./PROGRAMMES && chmod a+x ./PROGRAMMES/*
    fi

    # Run GONE
    bash ./script_GONE.sh ${pop} && \
        rm -rf TEMPORARY_FILES/

    # Prepare output
    awk -v pop=${pop} 'BEGIN {OFS="\\t"}; NR>2 { print pop, \$1, \$2}' Output_Ne_${pop} > ${pop}.ne
    """
}

process collectNe {
    label "small"
    publishDir "${params.outdir}/ne/results/", mode: 'copy', overwrite: true

    input:
    path pops

    output:
    path "GONE_NEVALS.ne"

    stub:
    """
    touch GONE_NEVALS.ne
    """

    script:
    """
    cat *.ne > GONE_NEVALS.ne
    """
}

