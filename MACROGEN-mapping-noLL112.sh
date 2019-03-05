#!/bin/bash

#####################################################################################
## Mapping MACROGEN FASTQ sequences (Illumina 1.9) to FELIS CATUS reference genome ##
#####################################################################################

# After testing with one single individual of Lynx lynx from the MACROGEN
# sequencing (LL112), this script will be used to go through the same steps
# (mapping, adding read groups and merging, marking duplicates and realigning) with
# the rest of the samples of the MACROGEN project.

# As explained in the 0.Mapping_pipeline.Rmd document, this mapping to the cat
# reference genome is necessary in order to generate demographic models through
# machine learning, a step which will be conducted by a collaborator (name ...),
# that will need high coverage sequencing data for at least 2 individuals per species.

# The Read Group Addition step has a starting code which will "extract" the exact run ID
# from the initial fastq list array. This does not apply for the MACROGEN samples, as
# their run ID only has one component. This part of the script is kept as it might apply
# to other samples I might use in the future.

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

# List of all Lynx lynx sample codes in MACROGEN project:
# (HERE I removed LL112 as the analysis for that sample were already run during testing)
MacroGenARRAY=($(ls /backup/grupolince/raw_data/MACROGEN/MACROGEN_trimmed/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1 | uniq | grep -v LL112))
# Path to cat reference genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
# No. of computer cores used. 20 = OK, >20 = ask people first!
THREADS=10
# Path to output files, were BAMS are generated:
OUT=/home/ebazzicalupo/CatRef_bams
# path to MACROGEN fastq files:
MacroGenPATH=/backup/grupolince/raw_data/MACROGEN/MACROGEN_trimmed/
# BARCODE MACROGEN, where the ID given by the sequencing company is converted to our sample ID
declare -A BARCODEID=(["LC1"]="c_lc_zz_0001" ["LL112"]="c_ll_vl_0112" ["LL146"]="c_ll_ya_0146" ["LL212"]="c_ll_cr_0212" ["LL90"]="c_ll_ki_0090" ["LR1"]="c_lr_zz_0001")

##########################
## Mapping with BWA MEM ##
##########################

for i in ${MacroGenARRAY[@]}
  do

    echo -e " *********** \n\n - Mapping ${i} - \n\n ***********"
    bwa mem $REF $MacroGenPATH/${i}_R1_trimmed.fastq.gz $MacroGenPATH/${i}_R2_trimmed.fastq.gz \
    -t $THREADS | samtools view -hbS -@ $THREADS - -o $OUT/${i}.cat_ref.bam
    echo " - Sorting ${i} -"
    samtools sort -@ $THREADS $OUT/${i}.cat_ref.bam -o $OUT/${i}.cat_ref.sorted.bam \
    && rm $OUT/${i}.cat_ref.bam

    echo -e " *********** \n\n - Adding READ Groups to ${i} and changing name to ${BARCODEID["${i}"]} - \n\n ***********"
    run=($(echo $i | cut -d"_" -f 1))  #Extracting run from i
    echo $run
    java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups \
    I=$OUT/${i}.cat_ref.sorted.bam \
    O=$OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg.bam \
    RGID=${i} RGLB=${BARCODEID["${i}"]}_lib \
    RGPL=Illumina RGPU=${run} RGSM=${BARCODEID["${i}"]} \
    VALIDATION_STRINGENCY=SILENT && rm $OUT/${i}.cat_ref.sorted.bam

    echo -e " *********** \n\n - Marking Duplicates of ${BARCODEID["${i}"]} and Re-Sorting - \n\n ***********"

    java -jar /opt/picard-tools/picard.jar MarkDuplicates \
    METRICS_FILE=$OUT/${i}_rmdup.txt \
    I=$OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg.bam \
    O=$OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup.bam \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
    rm $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg.bam
    samtools sort $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup.bam \
    -@ 10 -o $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted.bam
    rm $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup.bam
    samtools index $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted.bam

    echo -e " *********** \n\n - Realigning ${BARCODEID["${i}"]} - \n\n ***********"

    # RealignerTargetCreator
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -nt 10 -R $REF -I $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${BARCODEID["${i}"]}_realignertargetcreator.intervals
    # IndelRealigner
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $REF -targetIntervals $OUT/${BARCODEID["${i}"]}_realignertargetcreator.intervals \
    -I $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
    rm $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted.bam
    samtools flagstat $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    > $OUT/${BARCODEID["${i}"]}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.stats

done
