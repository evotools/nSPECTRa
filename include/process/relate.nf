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
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

    input:
    path vcf
    path tbi
    path ancfa
    path ancfai
    path maskfa
    path maskfai
    path poplabels
    path relate
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

    ${relate}/bin/RelateFileFormats --mode ConvertFromVcf -i ${contig}.RECODE \
        --haps ${contig}.INPUT.haps \
        --sample ${contig}.INPUT.sample && rm ${contig}.RECODE.vcf.gz

    ${relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
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
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

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
    label "large_long_smallmem"
    publishDir "${params.outdir}/relate/relate", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(contig), path(haps), path(sample), path(map), path(dist), path(annot)
    path relate

    output:
    path "relate_chr${contig}.anc"
    path "relate_chr${contig}.mut"
    val contig

    stub:
    """
    touch relate_chr${contig}.anc
    touch relate_chr${contig}.mut
    """

    script:
    // Define effective population size.
    def ne = params.neval ? "-N ${params.neval}" : "-N 30000"
    // Define executable
    def relate = task.cpus == 1 ? "${params.relate}/bin/Relate" : "${params.relate}/scripts/RelateParallel/RelateParallel.sh"
    // Define resources on core count.
    def cores = task.cpus == 1 ? "" : "--threads ${task.cpus}"
    def memory = params.relate_memory ? "--memory ${params.relate_memory}" : ""
    """
    ${relate} --mode All \
        --mode All \
        -m ${params.mutation_rate} \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        --annot ${annot} \
        --dist ${dist} \
        -o relate_chr${contig} ${ne} ${cores} ${memory}
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
    path relate

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
    ${relate}/bin/RelateMutationRate \
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
    label "medium_largemem"

    input:
    path contig_csv
    tuple path(anc), path(mut), path(dist) //, path(dist)
    path ancestral
    path ancestralfai
    path maskfa
    path maskfai
    path poplabels
    path relate

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
    ${relate}/bin/RelateMutationRate \
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
    path relate

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
    ${relate}/bin/RelateMutationRate \
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

process relate_mut_chr_finalise {
    label "medium_mem"
    publishDir "${params.outdir}/relate/mut_finalised_chr", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(contig), path(anc), path(mut), path(mutbin), path(oppbin)
    path contig_csv
    path poplabels
    path relate

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
    ${relate}/bin/RelateMutationRate \
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
    //tuple path(ancs), path(muts), path(mutbins), path(oppbins), path(rates)
    path contig_csv
    path poplabels
    path relate

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

    ${relate}/bin/RelateMutationRate \
                --mode Finalize \
                --poplabels ${poplabels} \
                --chr chroms.txt \
                --years_per_gen ${params.intergen_time} \
                -i relate_mut_ne \
                -o relate_mut_ne
    """
}


process relate_ne {
    publishDir "${params.outdir}/relate/ne", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    path ancs
    path muts
    path contig_csv
    path poplabels
    path relate

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
    ${relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i relate \
            --chr chroms.txt \
            --poplabels ${poplabels} \
            -m 1.25e-8 \
            --years_per_gen ${params.intergen_time} \
            --threads ${task.cpus} \
            -o relate_mut_ne 2> ne.err

    # Generate single-pop Ne values
    coal2ne relate_mut_ne.pairwise.coal > relate_mut_ne.pairwise.ne
    """
}


process relate_mut_chr_pop {
    label "medium_mem"

    input:
    tuple val(pop), path(popfile), val(idx), val(contig)
    tuple path(anc), path(mut), path(dist) //, path(dist)
    path ancestral
    path ancestralfai
    path maskfa
    path maskfai
    path poplabels
    path relate

    output:
    tuple val(pop), path("relate_mut_ne_${pop}_chr${contig}_mut.bin"), path("relate_mut_ne_${pop}_chr${contig}_opp.bin")

    stub:
    """
    touch relate_mut_${pop}_chr${contig}_mut.bin
    touch relate_mut_${pop}_chr${contig}_opp.bin
    """

    script:
    """
    # Extract single-chr fastas
    samtools faidx ${ancestral} ${contig} > anc.${contig}.fa
    samtools faidx ${maskfa} ${contig} > mask.${contig}.fa

    # Run relate mutation spectra
    ${relate}/bin/RelateMutationRate \
        --mode ForCategoryForPopForChromosome \
        --ancestor anc.${contig}.fa \
        --mask mask.${contig}.fa \
        --years_per_gen ${params.intergen_time} \
        --poplabels ${poplabels} \
        --pop_of_interest ${pop} \
        --mutcat ${baseDir}/assets/k3.mutcat \
        -i relate_mut_ne_chr${contig} \
        -o relate_mut_ne_${pop}_chr${contig} 2> ctx.${pop}.${contig}.err

    # Remove extra fasta files
    rm anc.${contig}.fa
    rm mask.${contig}.fa
    """
}

process relate_chr_pop_mut_finalise {
    label "medium_mem"
    publishDir "${params.outdir}/relate/mut_finalised_chr_pop", mode: "${params.publish_dir_mode}", overwrite: true

    input:
    tuple val(pop), path(muts), path(opps)
    path contig_csv
    path relate

    output:
    path "relate_mut_ne_${pop}*"
    //tuple path("relate_mut_ne_chr${contig}.anc*"), path("relate_mut_ne_chr${contig}.mut*"), path("relate_mut_ne_chr${contig}_mut.bin"), path("relate_mut_ne_chr${contig}_opp.bin"), path("relate_mut_ne_chr${contig}.rate")

    stub:
    """
    touch relate_mut_ne_${pop}_mut.bin
    touch relate_mut_ne_${pop}_opp.bin
    touch relate_mut_ne_${pop}.rate
    """

    script:
    """
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} | sort -nk1,1 > chroms.txt 

    ${relate}/bin/RelateMutationRate --mode SummarizeForGenomeForCategory --chr chroms.txt -o relate_mut_ne_${pop}
    ${relate}/bin/RelateMutationRate --mode FinalizeForCategory -i relate_mut_ne_${pop} -o relate_mut_ne_${pop}
    """
}


process relate_plot_pop {
    label "renv"
    publishDir "${params.outdir}/relate/plot", mode: "${params.publish_dir_mode}", overwrite: true
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

    input:
    path rates 
    path poplabels

    output:
    path "plot_mutrate.pdf"

    stub:
    """
    touch plot_mutrate.pdf
    """

    script:
    """
    relate_plot_pop ${baseDir}/assets/k3.mutcat ${params.intergen_time} ${poplabels}
    """
}

/*
Create modules to run relate Ne faster
*/

process ne_removetreeswithfewmutations {
    label "renv_many"
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

    input:
    path ancs
    path muts
    val contig

    output:
    tuple val(contig), path("relate_mut_ne_subset_chr${contig}.anc.gz"), path("relate_mut_ne_subset_chr${contig}.mut.gz")
    
    stub:
    """
    touch relate_mut_ne_subset_chr${contig}.anc.gz
    touch relate_mut_ne_subset_chr${contig}.mut.gz
    """

    script:
    """
    # Run relate mutation spectra    
    ${params.relate}/bin/RelateExtract \
        --mode RemoveTreesWithFewMutations \
        --threshold ${params.relate_mutation_threshold} \
        --anc relate_chr${contig}.anc \
        --mut relate_chr${contig}.mut \
        -o relate_mut_ne_subset_chr${contig} 2> relate_mut_ne_subset_chr${contig}.log
    gzip relate_mut_ne_subset_chr${contig}.anc
    gzip relate_mut_ne_subset_chr${contig}.mut
    """
}

process ne_ConstRate {
    label "renv_many"
    publishDir "${params.outdir}/relate/ne", mode: "${params.publish_dir_mode}", overwrite: true
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

    input:
    path ancs
    path muts

    output:
    path "relate_mut_ne.coal"
    path "relate_mut_ne.pdf"


    stub:
    """
    touch relate_mut_ne.coal
    touch relate_mut_ne.pdf
    """

    script:
    """
    # Get chromosome list
    awk 'BEGIN{FS=","};{print \$NF}' ${contig_csv} > chroms.txt 

    # Run relate mutation spectra    
    ${params.relate}/bin/RelateCoalescentRate \
        --mode GenerateConstCoalFile \
        --years_per_gen ${params.intergen_time} \
        -i 30000 \
        -o relate_mut_ne
    coal2ne relate_mut_ne.coal
    """
}

process ne_Iterate {
    label "renv_many"
    executor 'local'
    publishDir "${params.outdir}/relate/ne/it", mode: "${params.publish_dir_mode}", overwrite: true
    conda (params.enable_conda ? "${baseDir}/envs/r-environment.yml" : null)

    input:
    tuple val(contig), path(anc), path(mut)

    output:
    path "./OUT/*.anc"
    path "./OUT/*.mut"


    stub:
    """
    touch relate_mut_ne.coal
    touch relate_mut_ne.pdf
    """

    script:
    """
    runN=1
    maxN=${params.n_iterations}
    while [ "\$runN" -le "\$maxN" ]; do
        nextflow run ${params.baseDir}/include/iterate/relate_iterate.nf \
            --runN "\$runN" \
            --maxN "\$maxN" \
            --relate ${params.relate} \
            --work_path \$PWD
            --base_name relate_mut_ne_subset_chr${contig}
            --outdir \$PWD/OUT
            -profile ${workflow.profile}
        runN=\$((runN + 1))
    done
    """
}