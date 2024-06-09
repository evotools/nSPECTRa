#!/bin/bash
# Load modules first
module load singularity
module load anaconda
module load igmm/apps/vcftools
source activate DataPy3

# Sif image of hal
SIF=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/Singularity/hal.sif

# Input parameters
# while getopts ":rtoav" opt; do
#   case $opt in
#     r ) ref=${OPTARG} ;;
#     t ) tgt=${OPTARG} ;;
#     o ) outd=${OPTARG} ;;
#     a ) hal=${OPTARG} ;;
#     v ) vars=${OPTARG} ;;
#   esac
# done
REF=$1
TGT=$2
OUTD=$3
HAL=$4
VARS=$5
# Print
echo "REFERENCE: ${REF}"
echo "TARGET   : ${TGT}"
echo "OUTPUTDIR: ${OUTD}"
echo "HAL      : ${HAL}"
echo "VARIANTS : ${VARS}"

# Extract the alignments for the two genomes
echo "Make alignments file from hal"
singularity exec $SIF hal2maf \
    --refGenome $REF \
    --targetGenomes $TGT \
    --hdf5InMemory $HAL ${OUTD}/${REF}_${TGT}.maf

# Create ancestral genome for the targeted species
echo "Make ancestral genome"
python MAF2ANCFA.py -m ${OUTD}/${REF}_${TGT}.maf -o ${OUTD}/Ancestral_${REF} -t 4 -O fasta
sed "s/Anc.${REF}.//g" ${OUTD}/Ancestral_${REF}.fasta > ${OUTD}/Ancestral_${REF}.rename.fasta
samtools faidx ${OUTD}/Ancestral_${REF}.rename.fasta
awk 'BEGIN{OFS="\t"};{print $1,"0",$2}' ${OUTD}/Ancestral_${REF}.rename.fasta.fai > ${OUTD}/${REF}_regions.bed

# Make original annotation
## Generate bed of positions
echo "Make bed file of positions"
bcftools view -v snps -m2 -M2 ${VARS} | sed 's/_pilon//g' |\
        awk -v val=0 'BEGIN{OFS="\t"}; $1!~"#" {print $1,$2-1-val,$2+val,$1"_"$2"_"$4"_"$5}' | \
        sed 's/,/-/g' > ${OUTD}/${REF}_variants2lift_1nt.bed

## Lift position to ancestral
echo "Liftover positions"
singularity exec --bind $PWD:/mnt ${SIF} halLiftover \
        --hdf5InMemory $HAL \
        ${REF} ${OUTD}/${REF}_variants2lift_1nt.bed \
        ${TGT} ${OUTD}/${REF}_variantslifted_1nt.bed 2> ${OUTD}/MISSED_1nt.txt

# Get target fasta
echo "Make ${TGT} fasta"
singularity exec --bind $PWD:/mnt ${SIF} hal2fasta \
        --hdf5InMemory $HAL ${TGT} > ${OUTD}/${TGT}.fasta

## Extract base of interest, fix strandedness and create input for bcftools annotate
echo "Create annotation file for Variants"
bedtools getfasta -fi ${OUTD}/${TGT}.fasta -bed ${OUTD}/${REF}_variantslifted_1nt.bed -name -tab -fo ${OUTD}/Ancestral_${REF}.tab
./GetRightAllele.py ${OUTD}/Ancestral_${REF}.tab | \
    awk 'BEGIN{OFS="\t"}; NR>1{print $1,$2,$4}' | \
    bgzip -c > ${OUTD}/Ancestral_annotation_${REF}.txt.gz && \
    tabix -b 2 -e 2 -s 1 ${OUTD}/Ancestral_annotation_${REF}.txt.gz

# Filter individuals
if [[ $VARS == *.bcf ]]; then
        plink --cow --bcf $VARS --allow-extra-chr --rel-cutoff 0.025 --double-id --missing --out ${OUTD}/missingness_${REF}
else
        plink --cow --vcf $VARS --allow-extra-chr --rel-cutoff 0.025 --double-id --missing --out ${OUTD}/missingness_${REF}
fi
awk '{print $2}' ${OUTD}/missingness_${REF}.rel.ids > ${OUTD}/missingness_${REF}.txt

# Load mutyper environment
echo "Run mutyper"
source deactivate && source activate mutyperenv
echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > ${OUTD}/header_${REF}.txt
bcftools view -R ${OUTD}/${REF}_regions.bed -v snps -m2 -M2 -S ${OUTD}/missingness_${REF}.txt $VARS |
    sed 's/_pilon//g' |
    bcftools annotate -a ${OUTD}/Ancestral_annotation_${REF}.txt.gz -h ${OUTD}/header_${REF}.txt -c CHROM,POS,AA |
    vcftools --vcf - --minGQ 30 --recode --recode-INFO-all --stdout |
    mutyper variants ${OUTD}/Ancestral_${REF}.rename.fasta - |
    mutyper spectra - > ${OUTD}/mutationSpectra_${REF}.txt