#!/bin/bash

# On genomics-b I copied filter5 vcf of the old dataset from cesga
# (/mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/CatRef_vcfs/Oscar_v1):
screen -S depth_oscar_v3

cd /home/ebazzicalupo/CatRef_vcfs
script depth_oscar_v3.log

# Reference Genome:
REF=/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# AllIndividuals input VCF
INVCF=WholeGenome_cat_ref.filter5.vcf

# VCF Directory
OUTdir=/home/ebazzicalupo/CatRef_vcfs

# Depth per sample Table generated in R (depth_loop.R)
DPStable=depth_per_sample.csv

# Dataset array
datasetARRAY=($(cat ${DPStable} | cut -d',' -f1 | grep -v "pop" ))

for dataset in ${datasetARRAY[@]}
 do
  
  # (1) Extract VCF of individuals of dataset
  samplesARRAY=($(cat sample_dataset.table | grep "${dataset}" | cut -f1))

  /opt/gatk-4.1.0.0/gatk SelectVariants \
  -R ${REF} \
  -V ${INVCF} \
  $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
  -O ${OUTdir}/${dataset}.vcf

  # Minimum and Maximum depth values (in this case a min of 5x per sample = column 11, and a max of mean+1.5*sd = column 8)
  # see depth_loop.R for more detail on each column value
  max=$(grep -w ${dataset} ${DPStable} | cut -d',' -f8)
  min=$(grep -w ${dataset} ${DPStable} | cut -d',' -f11)
  echo "Maximum depth of ${dataset} is ${max}, Minimum depth is ${min}"

  # (2) extract the excessive low/high depth variant with BCFtools filter
  echo "extracting excessively low/high depth variants from ${dataset} VCF"
  bcftools filter -i "INFO/DP < ${min} || INFO/DP > ${max}" -Ov ${OUTdir}/${dataset}.vcf \
  > ${OUTdir}/${dataset}.applydepthfilter.vcf

done

cp ${INVCF} ${OUTdir}/WholeGenome_cat_ref.filter7-3.vcf
OUTVCF=${OUTdir}/WholeGenome_cat_ref.filter7-3.vcf

datasetARRAY=($(ls *.applydepthfilter.vcf | tr ' ' '\n' | cut -d'.' -f1 | sort -u))

# Remove depth variants from VCF
for dataset in ${datasetARRAY[@]}
 do
    echo "subtracting excessively low/high depth variants of ${dataset} from output VCF"
    bedtools subtract -a ${OUTVCF} \
    -b ${OUTdir}/${dataset}.applydepthfilter.vcf -header \
    > tmp && mv tmp ${OUTVCF}
done

# Get list of high depth and low depth variants separated
for dataset in ${datasetARRAY[@]}
 do
    echo "getting list excessively high depth variants of ${dataset}"
    max=$(grep -w ${dataset} ${DPStable} | cut -d',' -f8)
    min=$(grep -w ${dataset} ${DPStable} | cut -d',' -f11)
    echo "Maximum depth of ${dataset} is ${max}, Minimum depth is ${min}"

    bcftools filter -i "INFO/DP < ${min}" -Ov ${OUTdir}/${dataset}.applydepthfilter.vcf \
    > ${OUTdir}/${dataset}.applydepthfilter.lowdepth.vcf

    bcftools filter -i "INFO/DP > ${max}" -Ov ${OUTdir}/${dataset}.applydepthfilter.vcf \
    > ${OUTdir}/${dataset}.applydepthfilter.highdepth.vcf
done

# get ExcessHet table to calculate distribution in R
grep -vn "#" WholeGenome_cat_ref.filter7-3.vcf | grep -o -n -E 'ExcessHet=[[:digit:]]{1,3}\.?[[:digit:]]{0,8}' |
cut -d '=' -f2 > ExcessHet.table
# download to laptop
scp ebazzicalupo@genomics-b.ebd.csic.es:~/CatRef_vcfs/ExcessHet.table .
