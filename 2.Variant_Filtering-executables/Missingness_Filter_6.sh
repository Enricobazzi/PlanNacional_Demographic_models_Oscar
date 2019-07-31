#!/bin/bash
#SBATCH -t 2-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

#############################################################
START=$(date)
echo "General_Filter_1-2-3-4 SCRIPT for $1 starting : $START"
#############################################################

####################################
## Missingness Filtering launcher ##
####################################

# With this script I want to apply a filter based on data missingness to my
# dataset of 20 individuals, composed of 11 Lynx pardinus, 4 Lynx lynx,
# 3 Lynx rufus and 2 Lynx canadiensis. I will remove variants which are completely
# absent in all the individuals of a species, except for Lynx pardinus.
# Lynx pardinus, having more individuals, will be filtered for variants absent in
# at least 70% of the individuals (8 out of 11).

# To do so I will use the following softwares:

# BEDtools 2.28.0
# BCFtools 1.9

module load bedtools
module load bcftools

###################################
## VARIABLE and PATHS definition ##
###################################


######################
## Applying filters ##
######################

# List all sample names in a .namelist file:
ls $LUSTRE/test/CatRef_bams/*.bam | rev | cut -d'/' -f1 | rev | cut -d '_' -f1-4 | sort -u \
> $LUSTRE/test/CatRef_bams/all-samples.namelist

# List of species:
speciesARRAY=($(ls $LUSTRE/test/CatRef_bams/*.bam | rev | cut -d'/' -f1 | rev | cut -d '_' -f2 | sort -u))

# Create a copy of the VCF file with a new name that will be used to filter out
# excessively missing variants:
cp $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter5.vcf $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf

# for each species: create a .namelist file, divide the VCF by species, extract the
# excessive missingness variant, filter them from the new (copied) VCF file of all
# individuals

echo "Per-Species missingness variant filtering:" > $LUSTRE/test/missingness.variants.log

for species in ${speciesARRAY[@]}
  do

  echo "extracting $species names"
  grep $species $LUSTRE/test/CatRef_bams/all-samples.namelist > $LUSTRE/test/CatRef_bams/$species.namelist

  echo "filtering $species individuals from original VCF"
  bcftools view -S $LUSTRE/test/CatRef_bams/$species.namelist -Ov \
  $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter5.vcf \
  > $LUSTRE/test/CatRef_vcfs/$species_cat_ref.filter5.subset.vcf

  echo "extracting missing variants from $species VCF and filtering them out"
  if [ $species == lc ]
    then
    bcftools filter -i "N_MISSING >= 2" -Ov $species_cat_ref.filter5.subset.vcf \
    > $species_cat_ref.filter5.subset.missing.vcf

    LCmiss=$(grep -v "#" $species_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LC : $LCmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf \
    -b $species_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf

  elif [ $species == ll ]
    then
    bcftools filter -i "N_MISSING >= 4" -Ov $species_cat_ref.filter5.subset.vcf \
    > $species_cat_ref.filter5.subset.missing.vcf

    LLmiss=$(grep -v "#" $species_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LL : $LLmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf \
    -b $species_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf

  elif [ $species == lp ]
    then
    bcftools filter -i "N_MISSING >= 8" -Ov $species_cat_ref.filter5.subset.vcf \
    > $species_cat_ref.filter5.subset.missing.vcf

    LPmiss=$(grep -v "#" $species_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LP : $LPmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf \
    -b $species_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf

  elif [ $species == lr ]
    then
    bcftools filter -i "N_MISSING >= 3" -Ov $species_cat_ref.filter5.subset.vcf \
    > $species_cat_ref.filter5.subset.missing.vcf

    LRmiss=$(grep -v "#" $species_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LR : $LRmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf \
    -b $species_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/WholeGenome_cat_ref.filter6.vcf

  fi
done
