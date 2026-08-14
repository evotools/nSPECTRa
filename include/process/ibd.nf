// processes to calculate the regions in identity by descent across populations
process ibd {
    tag "ibd ${contig}"
    label "medium"

    input:
    tuple val(contig), path(vcf), path(tbi), path(map)
    path refinedibd

    output:
    tuple val(contig), path(vcf), path(tbi), path(map), path("IBD.${contig}.ibd.gz")

    script:
    """
    javamem=`python -c "import sys; maxmem=int(sys.argv[1]); print( maxmem - int(maxmem * .1) )" ${ task.memory.toGiga() }`
    java -Xmx\${javamem}G -jar ${refinedibd} ${params.refined_ibd_params} \
            gt=${vcf} \
            out=IBD.${contig} \
            chrom=${contig} \
            nthreads=${task.cpus}
    """

    stub:
    """
    touch ibd.${contig}.ibd.gz
    """
}

process make_map {

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    
    tuple val(chrom), path(vcf), path(tbi), path("${vcf.simpleName}.map")

    script:
    """
    bcftools query -f '%CHROM\t%ID\t%POS\n' ${vcf} | awk 'BEGIN{OFS="\t"}; {print \$1,\$2,\$3 * ${params.recombination_rate},\$3}'> ${vcf.simpleName}.map
    """

    stub:
    """
    touch ${vcf.simpleName}.map
    """
}

process merge_ibd {
    tag "ibd ${contig}"
    label "medium"

    input:
    tuple val(contig), path(vcf), path(tbi), path(map), path(ibd)
    path merge_ibd

    output:
    path "IBD.${contig}.merged.ibd.gz"

    script:
    """
    javamem=`python -c "import sys; maxmem=int(sys.argv[1]); print( maxmem - int(maxmem * .1) )" ${ task.memory.toGiga() }`
    zcat ${ibd} | java -Xmx\${javamem}G -jar ${merge_ibd} \
            ${vcf} \
            ${map} \
            ${params.merge_ibd_params} | bgzip -c > IBD.${contig}.merged.ibd.gz
    """

    stub:
    """
    touch IBD.${contig}.merged.ibd.gz
    """
}