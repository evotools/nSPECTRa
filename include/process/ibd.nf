// processes to calculate the regions in identity by descent across populations
process ibd {
    tag "ibd ${contig}"
    label "medium"

    input:
    path vcf
    path tbi
    path refinedibd
    val contig

    output:
    tuple val(contig), path("ibd.${contig}.ibd.gz")

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
    path vcf
    path tbi

    output:
    path "${vcf.simpleName}.map"

    script:
    """
    bcftools query -f '%CHROM\t%ID\t%POS\n' ${vcf} | awk 'BEGIN{OFS="\t"}; {print \$1,\$2,\$3/1000000,\$3}'> ${vcf.simpleName}.map
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
    tuple val(contig), path(ibd)
    path merge_ibd
    path vcf
    path tbi
    path map

    output:
    path "ibd.${contig}.ibd.gz"

    script:
    """
    javamem=`python -c "import sys; maxmem=int(sys.argv[1]); print( maxmem - int(maxmem * .1) )" ${ task.memory.toGiga() }`
    java -Xmx\${javamem}G -jar ${merge_ibd} \
            ${ibd} \
            ${vcf} \
            ${map} \
            ${params.merge_ibd_params} \
            IBD.${contig}.merge
    """

    stub:
    """
    touch ibd.${contig}.ibd.gz
    """
}