def helpMessage() {
  log.info '''
========================================================================
     .-') _   .-')     _ (`-.    ('-.             .-') _   _  .-')   
    ( OO ) ) ( OO ).  ( (OO  ) _(  OO)           (  OO) ) ( \\( -O )  
,--./ ,--,' (_)---\\_)_.`     \\(,------.   .-----./     '._ ,------.  
|   \\ |  |\\ /    _ |(__...--'' |  .---'  '  .--./|'--...__)|   /`. ' 
|    \\|  | )\\  :` `. |  /  | | |  |      |  |('-.'--.  .--'|  /  | | 
|  .     |/  '..`''.)|  |_.' |(|  '--.  /_) |OO  )  |  |   |  |_.' | 
|  |\\    |  .-._)   \\|  .___.' |  .--'  ||  |`-'|   |  |   |  .  '.' 
|  | \\   |  \\       /|  |      |  `---.(_'  '--'\\   |  |   |  |\\  \\  
`--'  `--'   `-----' `--'      `------'   `-----'   `--'   `--' '--' 
========================================================================
      '''


params {
    annotation = false
    help = false
    publish_dir_mode = 'copy'
}

  
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run evotools/nSpectr -profile test,standard

    Input parameters:
      --variants [file.vcf.gz]        Path to vcf[.gz]/bcf of input variants.
      --idx [file.vcf.gz.tbi]         Path to vcf[.gz].tbi/bcf.csi of input variants. 
      --hal [file.hal]                Hal file with ancestral genome.
      --reference [genomeID]          Genome id present into the hal archive
      --target [ancestralID]          Ancestral genome id present into the hal archive
      --pops_folder [PATH]            Path to folder with lists of samples in each population (relate workflow)
      --annotation [PATH]             Path to a bed file with the negative strands section (for mutyper) 
      --ancestral [PATH]              Provide pre-computed ancestral genome to use instead of the hal file
    
    Algorithms selection:
      --algorithm [name]              Name of the algorithm to use for the mutation spectra definition (Default: mutyper; options: mutyper, relate, sdm)
      --imputation [name]             Name of the algorithm to use to perform imputation (Default: shapeit4; options: beagle, shapeit4)
      --species [name]                Species ID for VEP (e.g. 'bos_taurus')

    Imputation refinement:
      --neval [int]                   Effective population size to use for the imputation (Default: null)
    
    Effective population size estimates
      --ne_subset [int]               Number of markers to consider for Ne estimation in GONE (Default: 50000)

    External software:
      --beagle [PATH]                 Path to beagle exe (if unspecified, will download version 21Apr21.304)
      --shapeit [PATH]                Path to shapeit4 exe
      --relate [PATH]                 Path to the relate software suite

    VEP custom annotation:
      --download_cache                Download VEP cache for the given species (Default: false; if no cache and GFF are specified, it will try to download anyway)
      --no_sift                       This option disable the "--sift b" option in VEP (Default: false)
      --vep_cache                     Run VEP using local cache into path (Default: false; mandatory)
      --custom_vep                    Run VEP with provided fasta and gff files (Default: false)
      --gff [FILE.gff.gz]             Specify annotation GFF file for VEP (must match chromosome code in the fasta and the hal, and be sorted and bgzipped)
      --mask [FILE.bed]               Specify a bed file with regions to masks in relate and mutyper (No default specified); 
                                      if not provided, extract soft-masked regions from the genome analysed.
    Identity by descent for refined-ibd:
      --compute_ibd                   Compute refined-ibd (Default: false)
      --refined_ibd_params            Parameters for refined-ibd run (Default: 'window=10 trim=0.015 length=0.2 lod=4')
      --merge_ibd_params              Parameters to use to combine consecutive IBD segments as number of discordant homozygotes allowed and 
                                      and maximum distance allowed in centiMorgan (Default: '1 0.6')
      --refinedibd                    Path to merge-ibd jar (No default, if unspecified will download version 17Jan20.102)
      --mergeibd                      Path to merge-ibd jar (No default, if unspecified will download version 17Jan20.102)
     
    Filtering parameters:
      --vcf_is_filtered               The VCF file provided is pre-filtered (skip filtering)
      --imputed                       The VCF file provided is pre-imputed (skip imputation)
      --extract [BED]                 Limit to variants in given regions (bed format) 
      --exclude [BED]                 Limit to variants not in given regions (bed format) 
      --filter                        Filter data after imputation
      --coding                        Limit analyses to exonic variants only
      --noncoding                     Limit analyses to intronic and intergenic variants only
      --min_pop_size                  Minimum population size to consider for relate analyses (Default: 5)
      --chr_list [file]               List of chromosomes to use for vcf analyses in format "N,chrID" (if not provided, extract from vcf file provided)

    Other parameters:
      --filter_vcf [PARAMS]           Filtering parameters to use when --filter is on (for bcftools e.g. "-i 'F_MISSING<0.1'")
      --k [values]                    Values of K to use for Mutyper, comma separated if multiple values are considered (Default: 3,5,7)
      --intergen_time [value]         Define intergeneration time for the species (Default is 28, as in relate)
      --mutation_rate                 Define mutation rate for the species (Default: 1.25e-8)
      --cactus_url                    Provide url for the desired cactus release
      --rename_hal_sequences [file]   Rename sequences extracted from the cactus archive using this file

    Constrained elements detection
      --gone                          Run GONE to calculate effective population size (Default: false)
      --phast                         Provide path to phast installation
      --hal4d                         Use hal4d to define the neutral regions
      --exon_bed                      Provide path to bed file with exonic regions for hal4d

    Other
      --max_memory                    Max memory available on the machine (Default: '128.GB')
      --max_cpus                      Max cpu available on the machine (Default: 16)
      --max_time                      Max runtime available on the machine (Default: '240.h')
      --scratch [true/false]          Run the workflow with the scratch mode on (e.g. will use temporary files in temporary folder of node in cluster; default: false)
      --outdir [file]                 The output directory where the results will be saved (Default: './results')
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)

    Workflow profile
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: standard, conda, docker, singularity, podman, eddie, eddie_conda, sge, uge
      """
}
