


process get_hal {
    label "compile"

    output:
    path "cactus/"

    stub:
    """
    mkdir cactus
    mkdir cactus/bin
    touch cactus/bin/hal2maf
    touch cactus/bin/hal2fasta
    touch cactus/bin/halPhyloPTrain.py
    touch cactus/bin/halPhyloPMP.py
    """

    script:
    """
    wget -O - ${params.cactus_url} > cactus.tar.gz
    tar xvfz cactus.tar.gz && mv cactus-*/ cactus/ && rm cactus.tar.gz
    """

// Download beagle
process get_beagle {
    label "compile"

    output:
    path "beagle.21Apr21.304.jar"

    stub:
    """
    touch beagle.21Apr21.304.jar
    """

    script:
    """
    wget https://faculty.washington.edu/browning/beagle/beagle.21Apr21.304.jar
    """

}

// Download refined-ibd
process get_ref_ibd {
    label "compile"

    output:
    path "refined-ibd.17Jan20.102.jar"

    stub:
    """
    touch refined-ibd.17Jan20.102.jar
    """

    script:
    """
    wget https://faculty.washington.edu/browning/refined-ibd/refined-ibd.17Jan20.102.jar
    """

}

process get_merge_ibd {
    label "compile"

    output:
    path "merge-ibd-segments.17Jan20.102.jar"

    stub:
    """
    touch merge-ibd-segments.17Jan20.102.jar
    """

    script:
    """
    wget https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar
    """

}


// Get GONE for Ne estimate

 process gone_get { 
    label "small"

    output: 
    path "GONE"
  
    script:
    """
    # Get GONE
    git clone https://github.com/esrud/GONE.git
    """
}

process get_vep_cache {
    label "vep"

    output: 
    path "vep_cache"

    stub:
    """
    mkdir vep_cache
    """

    script:
    """
    # Get VEP installer
    vepver=`vep | awk '\$1=="ensembl" {print \$NF}'`
    wget https://raw.githubusercontent.com/Ensembl/ensembl-vep/release/\${vepver}/INSTALL.pl

    # Install vep cache
    mkdir vep_cache
    perl INSTALL.pl -a c -n --NO_BIOPERL -l -c ./vep_cache -s ${params.species}
    """
    
}
