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

################################
## Variant Filtering launcher ##
################################

# With this script I want to apply all of the General Hard filters to my VCF dataset.
# These include:

# (1) Repetitive/Low mappability regions
# (2) Indels + Non-biallelic sites
# (3) Lynx genus wide exclusive substitutions from Felis catus
# (4) and (5) Hard quality filters, as GATK standard practices

# It will generate an output hard-filtered.vcf file and a log file reporting the number
# of variants removed and the number of variants left at each step.

# As these filters are independent of species and sequencing technology they can be
# applied to the complete VCF file directly. A detailed explenation for each
# filtering step is available in the 2.Variant_Filtering.md MarkDown file.

# This script will run on the finis terrae II server of CESGA using the following softwares:

# GATK 4.1.1.0
# BEDtools 2.28.0
# BCFtools 1.9

# The VCF file name (without the .vcf extension) must be defined while launching
# the script as such:

# ./General_Filter_1-2-3-4.sh <VCFfilename>

module load gatk/4.1.1.0
module load bedtools
module load bcftools

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Input VCF File:
INvcf=$LUSTRE/test/CatRef_vcfs/"$1".vcf.gz

# Output LOG File:
OUTlog="$1".filtering.log
echo "General Variant Filtering LOG" > $OUTlog

# BED File of Masked regions:
MASKbed=$LUSTRE/test/CatGenome_Masked_BEDs/Masked_Regions.bed

# Step 1 output VCF
ST1out="$1".filter1.vcf

# Step 2 output VCF
ST2out="$1".filter2.vcf

# Step 3 output VCF
ST3out="$1".filter3.vcf

# Step 4 output VCF
ST4out="$1".filter4.vcf

# Step 5 output VCF
ST5out="$1".filter5.vcf


############################################
## (1) Repetitive/Low mappability regions ##
############################################

# Print number of variants before filtering to log:
ST1start=$(zcat $INvcf | grep -v "#" | wc -l)
echo "Repetitive/Low mappability regions - initial number of variants : $ST1start" >> $OUTlog

# Apply the filter with BedTools subtract
zcat $INvcf | bedtools subtract -a - -b $MASKbed -header | uniq > $ST1out

# Print number of variants after filtering to log:
ST1end=$(grep -v "#" $ST1out | wc -l)
echo "Repetitive/Low mappability regions - final number of variants : $ST1end" >> $OUTlog

# Print number of variants filtered to log:
ST1filtered="$(echo "$ST1start - $ST1end" | bc)"
echo "Repetitive/Low mappability regions - number of variants filtered : $ST1filtered" >> $OUTlog


######################################
## (2) Indels + Non-biallelic sites ##
######################################

# Apply the filter with GATK SelectVariants
gatk SelectVariants \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -R $REF \
  -V $ST1out \
  -O $ST2out

# Print number of variants after filtering to log:
ST2end=$(cat $ST2out | grep -v "#" | wc -l)
echo "Indels + Non-biallelic sites - final number of variants : $ST2end" >> $OUTlog

# Print number of variants filtered to log:
ST2filtered="$(echo "$ST1end - $ST2end" | bc)"
echo "Repetitive/Low mappability regions - number of variants filtered : $ST2filtered" >> $OUTlog


##################################################################
## (3) Lynx genus wide exclusive substitutions from Felis catus ##
##################################################################

# Apply the filter with BCFtools view
bcftools view -e 'INFO/AF=1.00' $ST2out > $ST3out

# Print number of variants after filtering to log:
ST3end=$(cat $ST3out | grep -v "#" | wc -l)
echo "Lynx genus wide exclusive substitutions - final number of variants : $ST3end" >> $OUTlog

# Print number of variants filtered to log:
ST3filtered="$(echo "$ST2end - $ST3end" | bc)"
echo "Lynx genus wide exclusive substitutions - number of variants filtered : $ST3filtered" >> $OUTlog


######################################
## (4) and (5) Hard quality filters ##
######################################

# Filter all except for the RanksSums:
gatk SelectVariants \
  --selectExpressions "QUAL >= 30 && QD >= 2.0 && FS <= 60.0 && MQ >= 40.0" \
  -R $REF \
  -V $ST3out \
  -O $ST4out

# Print number of variants after filtering to log:
ST4end=$(cat $ST4out | grep -v "#" | wc -l)
echo "Hard quality except for the RanksSums - final number of variants : $ST4end" >> $OUTlog

# Print number of variants filtered to log:
ST4filtered="$(echo "$ST3end - $ST4end" | bc)"
echo "Hard quality except for the RanksSums - number of variants filtered : $ST4filtered" >> $OUTlog


# Filter RankSums with bcftools view:
bcftools view -e 'INFO/MQRankSum<-12.5 | INFO/ReadPosRankSum<-8.0' $ST4out > $ST5out

# Print number of variants after filtering to log:
ST5end=$(cat $ST5out | grep -v "#" | wc -l)
echo "RanksSums excluded by Hard quality - final number of variants : $ST5end" >> $OUTlog

# Print number of variants filtered to log:
ST5filtered="$(echo "$ST4end - $ST5end" | bc)"
echo "RanksSums excluded by Hard quality - number of variants filtered : $ST5filtered" >> $OUTlog


###########################################################
END=$(date)
echo "General_Filter_1-2-3-4 SCRIPT for $1 ended : $START"
###########################################################
