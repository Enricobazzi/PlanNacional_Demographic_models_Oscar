#!/bin/bash

####################################################################################
## Mapping Bobcat1 FASTQ sequences (Illumina 1.9) to FELIS CATUS reference genome ##
####################################################################################

# This script will be used to go through the same steps (mapping, adding read groups
# and merging, marking duplicates and realigning) with the samples of bobcat that Janeka sent.

# The fastqs are in PN2017 (/GRUPOS/grupolince/) directory
# As for the PGENOMA mapping, there will be two loops. One for
# mapping and adding read groups, and the other for merging and the rest of the steps.

# As explained in the 0.Mapping_pipeline.Rmd document, this mapping to the cat
# reference genome is necessary in order to generate demographic models through
# machine learning, a step which will be conducted by a collaborator (name Oscar Lao),
# that will need high coverage sequencing data for at least 2 individuals per species.

# The Read Group Addition step has a starting code which will "extract" the exact run ID
# from the initial fastq list array.

#######################################################
## REFERENCE GENOME dictionary creation and indexing ##
#######################################################

# has already been done with the following commands:
#
# bwa index /home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
#
# samtools faidx /home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
#
# java -jar /opt/picard-tools/picard.jar CreateSequenceDictionary R= /home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa O= /home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.dict

###################################
## VARIABLE and PATHS definition ##
###################################

# List of all Lynx rufus sample codes in PN2017/Bobcat1:
Bobcat1ARRAY=($(ls /GRUPOS/grupolince/PN2017/Bobcat1 | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3 | uniq))
# Path to cat reference genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
# No. of computer cores used. 20 = OK, >20 = ask people first!
THREADS=10
# Path to output files, were BAMS are generated:
OUT=/home/ebazzicalupo/CatRef_bams
# path to fastq files:
Bobcat1PATH=/GRUPOS/grupolince/PN2017/Bobcat1
# BARCODES, where the ID of the fastq file is converted to our sample ID:
declare -A BARCODEID=(["Bobcat1_S1_L001"]="c_lr_nm_0006" ["Bobcat1_S1_L002"]="c_lr_nm_0006" ["Bobcat1_S1_L003"]="c_lr_nm_0006" ["Bobcat1_S1_L004"]="c_lr_nm_0006")

##########################################
## Mapping Lynx rufus with BWA MEM ##
##########################################

for i in ${Bobcat1ARRAY[@]}
  do

    echo -e " ********************** \n\n - Mapping ${i} - \n\n **********************"
    bwa mem $REF $Bobcat1PATH/${i}_R1_001.fastq.gz $Bobcat1PATH/${i}_R2_001.fastq.gz \
    -t $THREADS | samtools view -hbS -@ $THREADS - -o $OUT/${i}.cat_ref.bam
    echo " - Sorting ${i} -"
    samtools sort -@ $THREADS $OUT/${i}.cat_ref.bam -o $OUT/${i}.cat_ref.sorted.bam \
    && rm $OUT/${i}.cat_ref.bam

    echo -e " ********************** \n\n - Adding READ Groups to ${i} and changing name to ${BARCODEID["${i}"]} - \n\n **********************"
    run=($(echo $i | cut -d"_" -f 1))  #Extracting run from i
    echo $run
    java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups \
    I=$OUT/${i}.cat_ref.sorted.bam \
    O=$OUT/${BARCODEID["${i}"]}_${i}_cat_ref_sorted_rg.bam \
    RGID=${i} RGLB=${BARCODEID["${i}"]}_lib \
    RGPL=Illumina RGPU=${run} RGSM=${BARCODEID["${i}"]} \
    VALIDATION_STRINGENCY=SILENT && rm $OUT/${i}.cat_ref.sorted.bam

done

LRSampleARRAY=($(for i in ${Bobcat1ARRAY[@]}; do echo ${BARCODEID["${i}"]}; done | sort -u))


for i in ${LRSampleARRAY[@]}
  do

    echo -e " ********************** \n\n - Merging all ${i} samples and Re-Sorting - \n\n **********************"

    ls $OUT/${i}_*_sorted_rg.bam  > $OUT/${i}.bam.list
    samtools merge  -@ $THREADS -r $OUT/${i}_merged.bam -b $OUT/"${i}".bam.list

    BAMARRAY=($(cat $OUT/${i}.bam.list))
    for k in ${BAMARRAY[@]}
      do
        echo "Removing $k"
        rm $k
    done

    echo " Re-Sorting ${i}"
    samtools sort  -@ $THREADS $OUT/${i}_merged.bam -o $OUT/${i}_cat_ref_sorted_rg.bam && rm $OUT/${i}_merged.bam;


    echo -e " ********************** \n\n - Marking Duplicates of ${i} and Re-Sorting - \n\n **********************"

    java -jar /opt/picard-tools/picard.jar MarkDuplicates \
    METRICS_FILE=$OUT/${i}_rmdup.txt \
    I=$OUT/${i}_cat_ref_sorted_rg.bam \
    O=$OUT/${i}_cat_ref_sorted_rg_rmdup.bam \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
    rm $OUT/${i}_cat_ref_sorted_rg.bam
    samtools sort $OUT/${i}_cat_ref_sorted_rg_rmdup.bam \
    -@ 10 -o $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam
    rm $OUT/${i}_cat_ref_sorted_rg_rmdup.bam
    samtools index $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam

    echo -e " ********************** \n\n - Realigning ${i} - \n\n **********************"

    # RealignerTargetCreator
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -nt 10 -R $REF -I $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${i}_realignertargetcreator.intervals
    # IndelRealigner
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $REF -targetIntervals $OUT/${i}_realignertargetcreator.intervals \
    -I $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
    rm $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam
    samtools flagstat $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    > $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.stats

done
