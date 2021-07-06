// make popfile for relate
process makepopfile {
    tag "popfile"
    label "small"

    input:
    path pop_folder

    output:
    path "popfile.poplabels"

    
    stub:
    """
    echo "sample population group sex" > popfile.poplabels
    echo "sample1 population2 group2 NA" >> popfile.poplabels
    echo "sample2 population1 group2 NA" >> popfile.poplabels
    echo "sample3 population1 group1 NA" >> popfile.poplabels
    """

    script:
    """
    makePopfile ${pop_folder} > popfile.poplabels
    """
}

// Get chromosome list
process chromosomeList {
    tag "chr"
    label "small"

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
    bcftools view ${ch_vcf} |\
        awk '\$1!~"#" {print \$1}' |\
        sort -T . |\
        uniq |\
        awk 'BEGIN{OFS=","};{print NR, \$0}' > sequences.csv 
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
    """
}

process make_relate_map {
    label "renv"

    input:
    tuple val(contig), path(haps), path(sample), path(dist), path(annot)

    output:
    tuple val(contig), path(haps), path(sample), path("${contig}.map"), path(dist), path(annot)
    
    stub:
    """
    touch ${contig}.map
    """

    script:
    """
    # Make genetic maps
    makeRelateMaps ${contig}
    """
}

// Relate main process
process relate {
    label "medium_mem"
    publishDir "${params.outdir}/relate/relate", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(contig), path(haps), path(sample), path(map), path(dist), path(annot)

    output:
    path "relate_chr${contig}.anc"
    path "relate_chr${contig}.mut"

    stub:
    """
    touch relate_chr${contig}.anc
    touch relate_chr${contig}.mut
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
        -o relate_chr${contig} 
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
        -o relate_chr${contig} 
    """
}

// Average mutation rate
process relate_avg_mut {
    label "medium_mem"
    publishDir "${params.outdir}/relate/average_mut", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple path(ancs), path(muts), path(dist)
    path contig_csv
    path poplabels

    output:
    path "relate_mut_gen_avg.rate"


    stub:
    """
    touch relate_mut_gen_avg.rate
    touch avg.gen.err
    touch avg.gen.out
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
                --poplabels ${poplabels} \
                -i relate_mut_ne \
                -o relate_mut_gen > avg.time.out 2> avg.time.err
    
    python -c "import sys; vals = [float(line.strip().split()[1]) for line in open(sys.argv[1]) if 'nan' not in line]; print(sum(vals)/ len(vals))" relate_mut_gen_avg.rate > mut_rate.txt
    """
}

process make_mask {
    label "medium"

    input:
    path genome
    path bedfile

    output:
    path "mask.fasta"

    stub:
    """
    touch mask.fasta
    """

    script:
    """
    bedtools maskfasta -fi ${genome} -bed ${bedfile} -fo tmp.fasta
    python -c 'import sys; [sys.stdout.write(line) if ">" in line else sys.stdout.write(line.upper().replace("A","P").replace("C","P").replace("G","P").replace("T","P").replace("-","N") ) for line in open(sys.argv[1])]' tmp.fasta > mask.fasta 
    rm tmp.fasta
    """
}

// Relate mutation 
process relate_mut {
    label "medium_mem"

    input:
    path contig_csv
    tuple path(anc), path(mut), path(dist) //, path(dist)
    path ancestral
    path ancestralfai
    path maskfa
    path maskfai
    path poplabels

    output:
    path "ANCMUT/relate_mut_ne*"
    
    stub:
    """
    mkdir ANCMUT
    cp relate_mut_ne_chr${contig}* ANCMUT/
    touch ANCMUT/relate_mut_ne_mut.bin
    touch ANCMUT/relate_mut_ne_opp.bin
    """

    script:
    """
    # Extract chrs
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} | sort -nk1,1 > chroms.txt 

    # Copy inputs locally
    mkdir ANCMUT
    cp relate_mut_ne* ANCMUT/

    # Run relate mutation spectra
    ${params.relate}/bin/RelateMutationRate \
                --mode WithContext \
                --ancestor ${ancestral} \
                --mask ${maskfa} \
                --chr chroms.txt \
                --years_per_gen ${params.intergen_time} \
                --poplabels ${poplabels} \
                -i ANCMUT/relate_mut_ne \
                -o ANCMUT/relate_mut_ne 2> ctx.err

    """
}

process relate_mut_chr {
    label "medium_mem"

    input:
    tuple val(idx), val(contig)
    tuple path(anc), path(mut), path(dist) //, path(dist)
    path ancestral
    path ancestralfai
    path maskfa
    path maskfai
    path poplabels

    output:
    tuple val(contig), path("ANCMUT/relate_mut_ne_chr${contig}.anc*"), path("ANCMUT/relate_mut_ne_chr${contig}.mut*"), path("ANCMUT/relate_mut_ne_chr${contig}_mut.bin"), path("ANCMUT/relate_mut_ne_chr${contig}_opp.bin")

    stub:
    """
    mkdir ANCMUT
    touch ANCMUT/relate_mut_chr${contig}.anc.gz
    touch ANCMUT/relate_mut_chr${contig}.mut.gz
    touch ANCMUT/relate_mut_chr${contig}_mut.bin
    touch ANCMUT/relate_mut_chr${contig}_opp.bin
    """

    script:
    """
    # Extract single-chr fastas
    samtools faidx ${ancestral} ${contig} > anc.${contig}.fa
    samtools faidx ${maskfa} ${contig} > mask.${contig}.fa

    # Copy anc/mut files locally
    mkdir ANCMUT
    cp relate_mut_ne_chr${contig}.* ANCMUT/

    # Run relate mutation spectra
    ${params.relate}/bin/RelateMutationRate \
                --mode WithContextForChromosome \
                --ancestor anc.${contig}.fa \
                --mask mask.${contig}.fa \
                --years_per_gen ${params.intergen_time} \
                --poplabels ${poplabels} \
                -i ANCMUT/relate_mut_ne_chr${contig} \
                -o ANCMUT/relate_mut_ne_chr${contig} 2> ctx.${contig}.err

    # Remove extra fasta files
    rm anc.${contig}.fa
    rm mask.${contig}.fa
    """
}

process relate_chr_mut_finalise {
    label "medium_mem"
    publishDir "${params.outdir}/relate/mut_finalised_chr", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(contig), path(anc), path(mut), path(mutbin), path(oppbin)
    path contig_csv
    path poplabels

    output:
    tuple path("FIN/relate_mut_ne_chr${contig}.anc*"), path("FIN/relate_mut_ne_chr${contig}.mut*"), path("FIN/relate_mut_ne_chr${contig}_mut.bin"), path("FIN/relate_mut_ne_chr${contig}_opp.bin"), path("FIN/relate_mut_ne_chr${contig}.rate")

    stub:
    """
    mkdir FIN
    touch FIN/relate_mut_ne_chr${contig}.anc.gz
    touch FIN/relate_mut_ne_chr${contig}.mut.gz
    touch FIN/relate_mut_ne_chr${contig}_mut.bin
    touch FIN/relate_mut_ne_chr${contig}_opp.bin
    touch FIN/relate_mut_ne_chr${contig}.rate
    """

    script:
    """
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} | sort -nk1,1 > chroms.txt 

    mkdir FIN
    cp relate_mut_ne_chr${contig}* FIN/
    ${params.relate}/bin/RelateMutationRate \
                --mode Finalize \
                --poplabels ${poplabels} \
                --years_per_gen ${params.intergen_time} \
                -i FIN/relate_mut_ne_chr${contig} \
                -o FIN/relate_mut_ne_chr${contig}
    """
}


process relate_mut_finalise {
    label "medium_mem"
    publishDir "${params.outdir}/relate/mut_finalised", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path inputs
    path contig_csv
    path poplabels

    output:
    path "relate_mut_ne.rate"
    path "relate_mut_ne_mut.bin"
    path "relate_mut_ne_opp.bin"

    stub:
    """
    touch relate_mut_ne.rate
    touch relate_mut_ne_mut.bin
    touch relate_mut_ne_opp.bin
    """

    script:
    """
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} | sort -nk1,1 > chroms.txt 

    ${params.relate}/bin/RelateMutationRate \
                --mode Finalize \
                --poplabels ${poplabels} \
                --chr chroms.txt \
                --years_per_gen ${params.intergen_time} \
                -i relate_mut_ne \
                -o relate_mut_ne
    """
}


process relate_ne {
    label "renv_multi"
    publishDir "${params.outdir}/relate/ne", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path ancs
    path muts
    path contig_csv
    path poplabels

    output:
    tuple path("relate_mut_ne_chr*.anc.gz"), path("relate_mut_ne_chr*.mut.gz"), path("relate_mut_ne_chr*.dist")
    tuple path("relate_mut_ne.coal"), path("relate_mut_ne.pairwise.coal"), path("relate_mut_ne.pairwise.bin")
    path "relate_mut_ne.pairwise.ne"
    path "relate_mut_ne.pdf"
    path "relate_mut_ne_avg.rate"


    stub:
    """
    for i in `awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv}`; do 
        touch relate_mut_ne_chr\${i}.anc.gz
        touch relate_mut_ne_chr\${i}.mut.gz
        touch relate_mut_ne_chr\${i}.dist
    done
    touch relate_mut_ne.coal
    touch relate_mut_ne.pairwise.coal
    touch relate_mut_ne.pairwise.bin
    touch relate_mut_ne_avg.rate
    """

    script:
    """
    # Get chromosome list
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} > chroms.txt 

    # Run relate mutation spectra    
    ${params.relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i relate \
            --chr chroms.txt \
            --poplabels ${poplabels} \
            -m 1.25e-8 \
            --years_per_gen ${params.intergen_time} \
            --threads ${task.cpus} \
            -o relate_mut_ne 2> ne.err

    # Generate single-pop Ne values
    coal2ne relate_mut_ne.pairwise.coal > relate_mut_ne.pairwise.ne

    #${params.relate}/bin/RelateCoalescentRate \
    #            --mode FinalizePopulationSize\
    #            --poplabels ${poplabels} \
    #            --chr chroms.txt \
    #            -i relate \
    #            -o relate 
    """
}

