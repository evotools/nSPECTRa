#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 * if Nsnps is set to 0, use all
 */


/*
 * Perform liftover from source to target using lastz (different species) 
 * or blat (same species).
 * Specify alignments options in the field alignerparam.
 * Outputs will be saved in the folder specified in outdir.
 * If an annotation is provided as bed/gff, liftover will lift it. 
 * If no annotation is provided, set to 'NO_FILE' or ''
 * Karyotype can be specified as ranges (e.g. 1-22), single different 
 * chromosomes can be added after comma (e.g. 1-22,X,Y,Mt).
 */
 
// Show help message
if (params.help) {
    include {helpMessage} from './include/process/help'
    helpMessage()
    exit 0
}


// Print run informations
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
log.info """\
Nextflow mutation Spectra v ${workflow.manifest.version}
========================================================================
variants        : $params.variants
idx             : $params.idx
output folder   : $params.outdir
hal             : $params.hal
reference       : $params.reference
target          : $params.target
procedure       : $params.algorithm
species         : $params.species
k               : $params.k
Ne subset       : $params.ne_subset
Intergen. time  : $params.intergen_time
Mut. rate       : $params.mutation_rate
Min. pop. size  : $params.min_pop_size
algorithm       : $params.algorithm
imputation sfw  : $params.imputation
filter          : $params.filter
coding          : $params.coding
noncoding       : $params.noncoding
annotation      : $params.annotation
pops_folder     : $params.pops_folder
pop labels file : $params.poplabels
chromosome list : $params.chr_list
cactus url      : $params.cactus_url 
phast           : $params.phast 
exons           : $params.exon_bed 
hal4d           : $params.hal4d 
relate          : $params.relate""" 
if (params.neval){
  log.info """Ne value        : $params.neval"""  
}
if (params.ancestral){
  log.info """ancestral       : $params.ancestral"""  
}
if (params.ref_fasta){
  log.info """refernece fasta : $params.ref_fasta"""  
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
checkPathParamList = [
    params.variants, params.idx, 
  ]
//    params.vcf, params.tbi ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/* Check mandatory params */
if (params.variants) { ch_var = file(params.variants) } else { exit 1, 'Vcf file not specified!' }
if (params.idx) { ch_var_idx = file(params.idx) } else { exit 1, 'TBI file not specified!' }
if (!params.hal) { exit 1, 'Hal file not specified and ancestral not specified!' }
if (!params.hal && !params.ref_fasta && !params.ancestral) { exit 1, 'Ancestral and reference genomes not specified!' }
if (params.hal4d && !params.exon_bed) { exit 1, 'Requested hal4d algorithm, but no bed with exons specified!' }
// algorithms = params.algorithm?.tokenize(',').flatten()

/*
 * Import sub-workflows
 */
include { PREPROCESS } from './include/workflow/vcfpreprocess' params(params)
include { ANCESTRAL } from './include/workflow/ancestral' params(params)
include { MUTYPER } from './include/workflow/mutyper' params(params)
include { RELATE } from './include/workflow/relate' params(params)
include { SDM } from './include/workflow/sdm' params(params)
include { GONE } from './include/workflow/gone' params(params)
include { IBD } from './include/workflow/ibd' params(params)
include { CONSTRAINED } from './include/workflow/phylop' params(params)
include { get_masks } from './include/process/prerun' 
include { get_hal } from './include/process/dependencies'

// Run workflows
workflow {
    // Generate the ancestral fasta
    get_hal()
    ANCESTRAL(get_hal.out)

    if (!params.ancestral_only){
      // Fetch mask files
      if (!params.mask){
        get_masks( ANCESTRAL.out[2] )
        ch_masks = get_masks.out
      } else {
        ch_masks = file(params.mask)
      }

      // Pre-process the variants of interest
      if (!params.vcf_is_filtered){
        PREPROCESS ( ch_var, ch_var_idx, ANCESTRAL.out[2], ANCESTRAL.out[3] )
        ch_var_new = PREPROCESS.out[0]
        ch_var_idx_new = PREPROCESS.out[1]

        // Get constrined elements and remove variants in them
        if ( params.phast || params.hal4d ){
          CONSTRAINED(get_hal.out, ch_var_new, ch_var_idx_new, PREPROCESS.out[2])
          ch_var_new = CONSTRAINED.out[0]
          ch_var_idx_new = CONSTRAINED.out[1]
        } 

        // Generate IBDs if requested
        if (params.compute_ibd){
          IBD()
        }
      } else {
        ch_var_new = ch_var
        ch_var_idx_new = ch_var_idx
      }

      // Run GONE to calculate Ne, if requested
      if (params.gone){
        GONE( ch_var_new, ch_var_idx_new )
      }

      // Run the actual mutation spectra
      if (params.algorithm =~ 'relate'){
          RELATE( ch_var_new, ch_var_idx_new, ANCESTRAL.out[0], ANCESTRAL.out[1], ANCESTRAL.out[2], ANCESTRAL.out[3], PREPROCESS.out[2], ch_masks )
      }
      if (params.algorithm =~ 'mutyper'){
          MUTYPER( ch_var_new, ch_var_idx_new, ANCESTRAL.out[0], ANCESTRAL.out[1], PREPROCESS.out[2], ch_masks )
          // MUTYPER( ch_var_new, ch_var_idx_new, ANCESTRAL.out[0], ANCESTRAL.out[1], PREPROCESS.out[2], ch_masks, GONE.out )
      } 
      if (params.algorithm =~ 'sdm'){
          SDM( ch_var_new, ch_var_idx_new, ANCESTRAL.out[0], ANCESTRAL.out[1], ANCESTRAL.out[2], ANCESTRAL.out[3], ch_masks, PREPROCESS.out[2] )
      }
    }

}