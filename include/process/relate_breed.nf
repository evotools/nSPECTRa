
process relate_format_legacy {
    label "renv"

    input:
    path vcf
    path tbi
    path ancfa
    path ancfai
    tuple val(samplename), path(samplelist), val(idx), val(contig)


    output:
    tuple val(samplename), val(contig), path("${contig}.${samplename}.RELATE.haps.gz"), path("${contig}.${samplename}.RELATE.sample.gz"), path("${contig}.map")
    
    stub:
    """
    touch ${contig}.${samplename}.RELATE.haps.gz
    touch ${contig}.${samplename}.RELATE.sample.gz
    touch ${contig}.map
    """

    script:
    """
    bcftools view -S ${samplelist} -r ${contig} --force-samples -O z ${vcf} > ${contig}.${samplename}.RECODE.vcf.gz

    ${params.relate}/bin/RelateFileFormats --mode ConvertFromVcf -i ${contig}.${samplename}.RECODE \
        --haps ${contig}.${samplename}.INPUT.haps \
        --sample ${contig}.${samplename}.INPUT.sample && rm ${contig}.${samplename}.RECODE.vcf.gz

    ${params.relate}/bin/RelateFileFormats \
                 --mode RemoveNonBiallelicSNPs \
                 --haps ${contig}.${samplename}.INPUT.haps \
                 -o ${contig}.${samplename}.RELATE_biallelic 

    ${params.relate}/bin/RelateFileFormats \
                 --mode FlipHapsUsingAncestor \
                 --haps ${contig}.${samplename}.RELATE_biallelic.haps \
                 --sample ${contig}.${samplename}.INPUT.sample \
                 --ancestor ${ancfa} \
                 -o ${contig}.${samplename}.RELATE_ancestral 

    cat ${contig}.${samplename}.INPUT.sample | gzip -c > ${contig}.${samplename}.RELATE.sample.gz
    cat ${contig}.${samplename}.RELATE_ancestral.haps | gzip -c > ${contig}.${samplename}.RELATE.haps.gz

    # Make genetic maps
    makeRelateMaps ${contig} ${samplename}

    #do some tidying
    rm ${contig}.${samplename}.INPUT.haps
    rm ${contig}.${samplename}.INPUT.sample
    rm ${contig}.${samplename}.RELATE_*
    """
}


// Relate format inputs by breed
process relate_format {
    label "renv"

    input:
    path vcf
    path tbi
    path ancfa
    path ancfai
    path maskfa
    path maskfai
    path poplabels
    tuple val(idx), val(contig)


    output:
    tuple val(contig), path("${contig}.RELATE.haps.gz"), path("${contig}.RELATE.sample.gz"), path("${contig}.RELATE.dist.gz"), path("${contig}.RELATE.annot")
    
    stub:
    """
    touch ${contig}.RELATE.haps.gz
    touch ${contig}.RELATE.sample.gz
    touch ${contig}.RELATE.annot
    touch ${contig}.RELATE.dist.gz
    touch ${contig}.map
    """

    script:
    """
    # Extract population of interest
    bcftools view -S ${samplelist} -r ${contig} --force-samples -O z ${vcf} > ${contig}.${samplename}.RECODE.vcf.gz

    # Extract single-chr fastas
    samtools faidx ${ancfa} ${contig} > anc.${contig}.fa
    samtools faidx ${maskfa} ${contig} > mask.${contig}.fa

    bcftools view -r ${contig} --force-samples -O z ${vcf} > ${contig}.RECODE.vcf.gz

    ${params.relate}/bin/RelateFileFormats --mode ConvertFromVcf -i ${contig}.RECODE \
        --haps ${contig}.INPUT.haps \
        --sample ${contig}.INPUT.sample && rm ${contig}.RECODE.vcf.gz

    ${params.relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
                    --haps ${contig}.INPUT.haps \
                    --sample ${contig}.INPUT.sample \
                    --ancestor anc.${contig}.fa \
                    --mask mask.${contig}.fa \
                    --poplabels ${poplabels} \
                    -o ${contig}.RELATE 
    
    # Make genetic maps
    makeRelateMaps ${contig}
    """
}


// Relate format inputs by breed
process relate_format_breed {
    label "renv"

    input:
    path vcf
    path tbi
    path ancfa
    path ancfai
    path maskfa
    path maskfai
    path poplabels
    tuple val(samplename), path(samplelist), val(idx), val(contig)


    output:
    tuple val(samplename), val(contig), path("${contig}.${samplename}.RELATE.haps.gz"), path("${contig}.${samplename}.RELATE.sample.gz"), path("${contig}.map"), path("${contig}.${samplename}.RELATE.dist.gz"), path("${contig}.${samplename}.RELATE.annot")
    
    stub:
    """
    touch ${contig}.${samplename}.RELATE.haps.gz
    touch ${contig}.${samplename}.RELATE.sample.gz
    touch ${contig}.${samplename}.RELATE.annot
    touch ${contig}.${samplename}.RELATE.dist.gz
    touch ${contig}.map
    """

    script:
    """
    # Extract single-chr fastas
    samtools faidx ${ancfa} ${contig} > anc.${contig}.fa
    samtools faidx ${maskfa} ${contig} > mask.${contig}.fa

    bcftools view -S ${samplelist} -r ${contig} --force-samples -O z ${vcf} > ${contig}.${samplename}.RECODE.vcf.gz

    ${params.relate}/bin/RelateFileFormats --mode ConvertFromVcf -i ${contig}.${samplename}.RECODE \
        --haps ${contig}.${samplename}.INPUT.haps \
        --sample ${contig}.${samplename}.INPUT.sample && rm ${contig}.${samplename}.RECODE.vcf.gz

    ${params.relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
                    --haps ${contig}.${samplename}.INPUT.haps \
                    --sample ${contig}.${samplename}.INPUT.sample \
                    --ancestor anc.${contig}.fa \
                    --mask mask.${contig}.fa \
                    --poplabels ${poplabels} \
                    -o ${contig}.${samplename}.RELATE 

    # Make genetic maps
    makeRelateMaps ${contig} ${samplename}

    #do some tidying
    rm ${contig}.${samplename}.INPUT.haps
    rm ${contig}.${samplename}.INPUT.sample
    rm anc.${contig}.fa
    rm mask.${contig}.fa
    """
}


// Relate main process by breed
process relate_br {
    label "medium_mem"
    publishDir "${params.outdir}/relate/relate/${samplename}", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(samplename), val(contig), path(haps), path(sample), path(map), path(dist), path(annot)

    output:
    tuple val(samplename), val(contig), path("relate_${samplename}_chr${contig}.anc"), path("relate_${samplename}_chr${contig}.mut")

    stub:
    """
    touch relate_${samplename}_chr${contig}.anc
    touch relate_${samplename}_chr${contig}.mut
    """

    script:
    if (task.cpus == 1)
    """
    ${params.relate}/bin/Relate --mode All \
        --mode All \
        -m ${params.mutation_rate} \
        -N ${params.neval} \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        --annot ${annot} \
        --dist ${dist} \
        --memory ${task.memory} \
        -o relate_${samplename}_chr${contig} 
    """
    else 
    """
    ${params.relate}/scripts/RelateParallel/RelateParallel.sh \
        --mode All \
        -m ${params.mutation_rate} \
        -N ${params.neval} \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        --annot ${annot} \
        --dist ${dist} \
        --threads ${task.cpus} \
        --memory ${task.cpus * task.memory} \
        -o relate_${samplename}_chr${contig} 
    """
}


// Average mutation rate by breed
process relate_avg_mut_br {
    label "medium_mem"
    publishDir "${params.outdir}/relate/average_mut/${samplename}", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(samplename), val(contigs), path(ancs), path(muts)
    path contig_csv

    output:
    path "relate_mut_${samplename}_gen_avg.rate"


    stub:
    """
    touch relate_mut_${samplename}_gen_avg.rate
    touch avg.gen.${samplename}.err
    """

    script:
    """
    # Get chromosome list
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} > chroms.txt 

    # Run average mutation rate
    ${params.relate}/bin/RelateMutationRate \
                --mode Avg \
                --chr chroms.txt \
                --years_per_gen ${params.intergen_time} \
                -i relate_${samplename} \
                -o relate_mut_${samplename}_gen > avg.time.${samplename}.out 2> avg.time.${samplename}.err
    """
}


// Relate mutation by breed
process relate_mut_br {
    label "medium_mem"

    input:
    tuple val(samplename), val(contig), path(anc), path(mut)//, path(dist)
    path ancestral
    path ancestralfai
    path maskfa
    path maskfai

    output:
    tuple val(samplename), val(contig), path(anc), path(mut), path("./relate_mut_${samplename}_chr${contig}_mut.bin"), path("./relate_mut_${samplename}_chr${contig}_opp.bin")

    stub:
    """
    touch relate_mut_${samplename}_chr${contig}_mut.bin
    touch relate_mut_${samplename}_chr${contig}_opp.bin
    """

    script:
    """
    # Extract single-chr fastas
    samtools faidx ${ancestral} ${contig} > anc.${contig}.fa
    samtools faidx ${maskfa} ${contig} > mask.${contig}.fa

    # Run relate mutation spectra
    ${params.relate}/bin/RelateMutationRate \
                 --mode WithContextForChromosome \
                 --ancestor anc.${contig}.fa \
                 --mask mask.${contig}.fa \
                 --years_per_gen ${params.intergen_time} \
                 -i relate_${samplename}_chr${contig} \
                 -o relate_mut_${samplename}_chr${contig} 2> ctx.${samplename}.${contig}.err

    # Remove extra fasta files
    rm anc.${contig}.fa
    rm mask.${contig}.fa
    """
}

