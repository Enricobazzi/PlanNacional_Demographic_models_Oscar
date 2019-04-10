#!/bin/bash

##########################################################################################
## Mapping American Lynx FASTQ sequences (Illumina 1.9) to FELIS CATUS reference genome ##
##########################################################################################

# After the MACROGEN and PGENOMA mapping, this script will be used to go through
# the same steps (mapping, adding read groups and merging, marking duplicates
# and realigning) with the samples of the American lynxes, that Murphy sent.

# As the fastqs are in two different folders the loops will be run separately for,
# each folder. The folders are in PN2017 (/GRUPOS/grupolince/)
# As for the PGENOMA mapping, there will be two loops (for each folder). One for
# mapping and adding read groups, and the other for merging and the rest of the steps.

# As explained in the 0.Mapping_pipeline.Rmd document, this mapping to the cat
# reference genome is necessary in order to generate demographic models through
# machine learning, a step which will be conducted by a collaborator (name ...),
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

# List of all Lynx canadiens 3 (Murphy) sample codes in PN2017:
LCA_3ARRAY=($(ls /GRUPOS/grupolince/PN2017/LCA_3/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2 | uniq))
# List of all Lynx rufus 30 (Murphy) sample codes in PN2017:
LRU_30ARRAY=($(ls /GRUPOS/grupolince/PN2017/LRU_30/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2 | uniq))
# Path to cat reference genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
# No. of computer cores used. 20 = OK, >20 = ask people first!
THREADS=10
# Path to output files, were BAMS are generated:
OUT=/home/ebazzicalupo/CatRef_bams
# path to fastq files:
LCA_3PATH=/GRUPOS/grupolince/PN2017/LCA_3
LRU_30PATH=/GRUPOS/grupolince/PN2017/LRU_30
# BARCODES, where the ID of the fastq file is converted to our sample ID:
declare -A BARCODEID=(["LCA-3_1"]="c_lc_zz_0003" ["LCA-3_2"]="c_lc_zz_0003" ["LCA-3_3"]="c_lc_zz_0003" ["LCA-3_4"]="c_lc_zz_0003" ["LCA-3_5"]="c_lc_zz_0003" ["LCA-3_6"]="c_lc_zz_0003" ["LCA-3_7"]="c_lc_zz_0003" ["LCA-3_8"]="c_lc_zz_0003" ["LRU-30_1"]="c_lr_fl_0005" ["LRU-30_2"]="c_lr_fl_0005" ["LRU-30_3"]="c_lr_fl_0005" ["LRU-30_4"]="c_lr_fl_0005" ["LRU-30_5"]="c_lr_fl_0005" ["LRU-30_6"]="c_lr_fl_0005" ["LRU-30_7"]="c_lr_fl_0005" ["LRU-30_8"]="c_lr_fl_0005")

##########################################
## Mapping Lynx canadensis with BWA MEM ##
##########################################

for i in ${LCA_3ARRAY[@]}
  do

    echo -e " ********************** \n\n - Mapping ${i} - \n\n **********************"
    bwa mem $REF $LCA_3PATH/${i}_AGTTCC_R1.fastq.gz $LCA_3PATH/${i}_AGTTCC_R2.fastq.gz \
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

LCSampleARRAY=($(ls $OUT/*c_lc_zz_0003*.bam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | sort -u))

for i in ${LCSampleARRAY[@]}
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

##########################################
## Mapping Lynx rufus with BWA MEM ##
##########################################

for i in ${LRU_30ARRAY[@]}
  do

    echo -e " ********************** \n\n - Mapping ${i} - \n\n **********************"
    bwa mem $REF $LRU_30PATH/${i}_GGTAGC_R1.fastq.gz $LRU_30PATH/${i}_GGTAGC_R2.fastq.gz \
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

LRSampleARRAY=($(ls $OUT/*c_lr_fl_0005*.bam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | sort -u))

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
