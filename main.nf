#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Import sub-workflows
 */
include {helpMessage} from './include/process/help'
include { PREPROCESS } from './include/workflow/vcfpreprocess'
include { ANCESTRAL } from './include/workflow/ancestral'
include { MUTYPER } from './include/workflow/mutyper'
include { RELATE } from './include/workflow/relate'
include { SDM } from './include/workflow/sdm'
include { GONE } from './include/workflow/gone'
include { IBD } from './include/workflow/ibd'
include { CONSTRAINED } from './include/workflow/phylop'
include { get_masks } from './include/process/prerun'
include { chromosomeList } from './include/process/prerun'

// Run workflows
workflow {
  // Show help message
  if (params.help) {
      helpMessage()
      exit 0
  }

  // Print run informations
  log.info '''\
  =================================================================================
      .-') _   .-')     _ (`-.    ('-.             .-') _   _  .-')     ('-.     
      ( OO ) ) ( OO ).  ( (OO  ) _(  OO)           (  OO) ) ( \\( -O )   ( OO ).-. 
  ,--./ ,--,' (_)---\\_)_.`     \\(,------.   .-----./     '._ ,------.   / . --. / 
  |   \\ |  |\\ /    _ |(__...--'' |  .---'  '  .--./|'--...__)|   /`. '  | \\-.  \\  
  |    \\|  | )\\  :` `. |  /  | | |  |      |  |('-.'--.  .--'|  /  | |.-'-'  |  | 
  |  .     |/  '..`''.)|  |_.' |(|  '--.  /_) |OO  )  |  |   |  |_.' | \\| |_.'  | 
  |  |\\    |  .-._)   \\|  .___.' |  .--'  ||  |`-'|   |  |   |  .  '.'  |  .-.  | 
  |  | \\   |  \\       /|  |      |  `---.(_'  '--'\\   |  |   |  |\\  \\   |  | |  | 
  `--'  `--'   `-----' `--'      `------'   `-----'   `--'   `--' '--'  `--' `--' 
  =================================================================================
  '''

  log.info """\
  Nextflow Mutation Spectra v${workflow.manifest.version}
  =================================================================================
  variants        : $params.variants
  idx             : $params.idx
  output folder   : $params.outdir
  hal             : $params.hal
  reference       : $params.reference
  target          : $params.target
  mutyper         : $params.mutyper
  sdm             : $params.sdm
  relate          : $params.relate
  species         : $params.species
  k               : $params.k
  Ne subset       : $params.ne_subset
  Intergen. time  : $params.intergen_time
  Mut. rate       : $params.mutation_rate
  Min. pop. size  : $params.min_pop_size
  imputation sfw  : $params.imputation
  coding          : $params.coding
  noncoding       : $params.noncoding
  annotation      : $params.annotation
  pops_folder     : $params.pops_folder
  pop labels file : $params.poplabels
  chromosome list : $params.chr_list
  cactus url      : $params.cactus_url 
  exons           : $params.exon_bed 
  constrained     : $params.constrained 
  relate dir      : $params.relate""" 
  if (params.neval){
    log.info """Ne value        : $params.neval"""  
  }
  if (params.ancestral_fna){
    log.info """ancestral       : $params.ancestral_fna"""  
  }
  if (params.reference_fna){
    log.info """reference fasta : $params.reference_fna"""  
  }
  if (params.mask){
    log.info """mask            : $params.mask"""  
  }
  if (params.beagle){
    log.info """beagle path     : $params.beagle"""  
  }
  if (params.shapeit){
    log.info """shapeit path    : $params.shapeit"""  
  }
  if (params.custom_vep){
    log.info """VEP gff annot.  : $params.gff"""  
  } else {
    log.info """VEP cache       : $params.vep_cache"""  
  }
  if (params.ancestral_only){
    log.info """ 
    
    Running generation of ancestral genome only
    
    """  
  }

  // Check parameters
  Channel.from([params.variants, params.idx])
  | map{
    if (!file(it).exists()) {
      throw new Exception("Missing parameter ${it}")
    }
  }

  /* Check mandatory params */
  if (params.variants && !params.ancestral_only) { ch_var = Channel.fromPath(params.variants) } else { exit 1, 'Vcf file not specified!' }
  if (params.idx && !params.ancestral_only) { ch_var_idx = Channel.fromPath(params.idx) } else { exit 1, 'TBI file not specified!' }
  if (!params.hal) { exit 1, 'Hal file not specified and ancestral not specified!' }
  if (!params.hal && !params.reference_fna && !params.ancestral_fna) { exit 1, 'Ancestral and reference genomes not specified!' }
  if (params.constrained && !params.exon_bed) { exit 1, 'Requested hal4d algorithm, but no bed with exons specified!' }
  if (!params.species){ throw new Exception("Parameter --species is required for file naming.") }

  // Generate the ancestral fasta
  ANCESTRAL()
  reference_fna = ANCESTRAL.out.ref_fna
  reference_fai = ANCESTRAL.out.ref_fai
  ancestral_fna = ANCESTRAL.out.anc_fna
  ancestral_fai = ANCESTRAL.out.anc_fai

  if (!params.ancestral_only){
    // Fetch mask files
    if (!params.mask){
      get_masks( reference_fna )
      ch_masks = get_masks.out
    } else {
      ch_masks = file(params.mask)
    }

    // Pre-process the variants of interest
    if (!params.vcf_is_filtered){
      PREPROCESS ( ch_var, ch_var_idx, reference_fna, reference_fai, ancestral_fna, ancestral_fai, ch_masks )
      ch_var_new = PREPROCESS.out.vcf
      ch_var_idx_new = PREPROCESS.out.tbi
      ch_chr_lists = PREPROCESS.out.chroms
      vcf_by_chr = PREPROCESS.out.vcf_by_chr
      vcf_chunks_ch = PREPROCESS.out.chunks_ch

      // Get constrined elements and remove variants in them
      if ( params.constrained ){
        CONSTRAINED(ch_var_new, ch_var_idx_new, ch_chr_lists)
        ch_var_new = CONSTRAINED.out[0]
        ch_var_idx_new = CONSTRAINED.out[1]
      } 

      // Generate IBDs if requested
      if (params.compute_ibd){
        IBD(ch_var_new, ch_var_idx_new, ch_chr_lists)
      }
    } else {
      ch_var_new = ch_var
      ch_var_idx_new = ch_var_idx
      chromosomeList( ch_var, ch_var_idx )
      ch_chr_lists = chromosomeList.out
      | combine(
        ch_var_new
        | combine( ch_var_idx )
      )
    }

    // Run GONE to calculate Ne, if requested
    if (params.gone){
      GONE( ch_var_new, ch_var_idx_new )
    }

    // Run the actual mutation spectra
    if (params.relate){
        RELATE( ch_var_new, ch_var_idx_new, ancestral_fna, ancestral_fai, ch_chr_lists, ch_masks )
    }
    if (params.mutyper){
        MUTYPER( vcf_by_chr, ch_var_idx_new, ancestral_fna, ancestral_fai, ch_chr_lists, ch_masks, vcf_chunks_ch )
    } 
    if (params.sdm){
        SDM( vcf_by_chr, reference_fna, reference_fai, ch_masks, ch_chr_lists )
    }
  }
}