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
echo "Missingness_Filter_6 SCRIPT for $1 starting : $START"
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

# # The VCF file name (without the extensions) must be defined while launching
# the script as such:

# ./Missingness_Filter_6.sh <VCFfilename>

#####################################
## Applying filters - Preparations ##
#####################################

# List all sample names in a .namelist file:
ls $LUSTRE/test/CatRef_bams/*.bam | rev | cut -d'/' -f1 | rev | cut -d '_' -f1-4 | sort -u \
> $LUSTRE/test/CatRef_bams/all-samples.namelist

# List species in an array (for loop):
speciesARRAY=($(ls $LUSTRE/test/CatRef_bams/*.bam | rev | cut -d'/' -f1 | rev | cut -d '_' -f2 | sort -u))

# Create a copy of the VCF file with a new name that will be used to filter out
# excessively missing variants:
cp $LUSTRE/test/CatRef_vcfs/"$1".filter5.vcf $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

#######################################
## Applying filters - The great LOOP ##
#######################################

# For each species:
# (1) create a .namelist file (with the names of all samples of that species);
# (2) use the namelist file to divide the VCF by species;
# (3) extract the excessive missingness variant with BCFtools filter, F_MISSING indicates
#       the proportion of missing data, the thresholds are explained above (lines 19-24);
# (4) filter the excessively missing variants from the new (lines 50-52) VCF file of all samples
#       with BEDtools subtract.

# Have a log file with filtered variants counts:
echo "Per-Species missingness variant filtering:" > $LUSTRE/test/missingness.variants.log

for species in ${speciesARRAY[@]}
  do

# (1) create a .namelist file (with the names of all samples of that species)
  echo "extracting $species names"
  grep $species $LUSTRE/test/CatRef_bams/all-samples.namelist > $LUSTRE/test/CatRef_bams/$species.namelist

# (2) use the namelist file to divide the VCF by species
  echo "filtering $species individuals from original VCF"
  bcftools view -S $LUSTRE/test/CatRef_bams/"$species".namelist -Ov \
  $LUSTRE/test/CatRef_vcfs/"$1".filter5.vcf \
  > $LUSTRE/test/CatRef_vcfs/"$species"_cat_ref.filter5.subset.vcf

# (3) extract the excessive missingness variant with BCFtools filter and
# (4) filter the excessively missing variants from the new VCF file of all samples
  echo "extracting missing variants from $species VCF and filtering them out"
  if [ $species == lc ]
    then
    bcftools filter -i "F_MISSING = 1" -Ov "$species"_cat_ref.filter5.subset.vcf \
    > "$species"_cat_ref.filter5.subset.missing.vcf

    LCmiss=$(grep -v "#" "$species"_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LC : $LCmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf \
    -b "$species"_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

  elif [ $species == ll ]
    then
    bcftools filter -i "F_MISSING = 1" -Ov "$species"_cat_ref.filter5.subset.vcf \
    > "$species"_cat_ref.filter5.subset.missing.vcf

    LLmiss=$(grep -v "#" "$species"_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LL : $LLmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf \
    -b "$species"_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

  elif [ $species == lp ]
    then
    bcftools filter -i "F_MISSING = 0.7" -Ov "$species"_cat_ref.filter5.subset.vcf \
    > "$species"_cat_ref.filter5.subset.missing.vcf

    LPmiss=$(grep -v "#" "$species"_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LP : $LPmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf \
    -b "$species"_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

  elif [ $species == lr ]
    then
    bcftools filter -i "F_MISSING = 1" -Ov "$species"_cat_ref.filter5.subset.vcf \
    > "$species"_cat_ref.filter5.subset.missing.vcf

    LRmiss=$(grep -v "#" "$species"_cat_ref.filter5.subset.missing.vcf | wc -l)
    echo "Variants filtered for LR : $LRmiss" >> $LUSTRE/test/missingness.variants.log

    bedtools subtract -a $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf \
    -b "$species"_cat_ref.filter5.subset.missing.vcf -header \
    > tmp && mv tmp $LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

  fi
done

###########################################################
END=$(date)
echo "Missingness_Filter_6 SCRIPT for $1 ended : $END"
###########################################################
