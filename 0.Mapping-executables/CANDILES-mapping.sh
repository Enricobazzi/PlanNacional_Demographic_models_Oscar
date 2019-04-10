#!/bin/bash

###############################################################################
## Mapping CANDILES sequences (Illumina 1.5) to FELIS CATUS reference genome ##
###############################################################################

# This script will be used to go through the same steps (mapping, adding read groups and merging,
# marking duplicates and realigning) with the remaining sample CANDILES (Lypa23).

# The fastq files were transfered to the genomics-b.ebd.csic.es server from
# from the cesga server in a folder called Enrico_moveLypa with the following path:
# /home/ebazzicalupo/fastqs/Enrico_moveLypa

# tar -zcvf Enrico_moveLypa.tar.gz Enrico_moveLypa
# scp Enrico_moveLypa.tar.gz ebazzicalupo@genomics-b.ebd.csic.es:~
# tar -zxvf Enrico_moveLypa.tar.gz

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

# List of all Lynx pardinus sample codes in Lypa23:
CandilesARRAY=($(ls /home/ebazzicalupo/fastqs/Enrico_moveLypa/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | uniq))
# Path to cat reference genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
# No. of computer cores used. 20 = OK, >20 = ask people first!
THREADS=10
# Path to output files, were BAMS are generated:
OUT=/home/ebazzicalupo/CatRef_bams
# path to Candiles fastq files:
CandilesPATH=/home/ebazzicalupo/fastqs/Enrico_moveLypa
# BARCODE CANDILES copied from Maria:
declare -A BARCODEID=(["6220RAAXX_lane3_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane4_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane6_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane7_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane8_sequence_0"]="c_lp_sm_0221" ["62AHEAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["621CYAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane5_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane2_sequence_0"]="c_lp_sm_0221")

##########################
## Mapping with BWA MEM ##
##########################

for i in ${CandilesARRAY[@]}
  do

    echo -e " ********************** \n\n - Mapping ${i} - \n\n **********************"
    bwa mem $REF $CandilesPATH/${i}_1.fastq.gz $CandilesPATH/${i}_2.fastq.gz \
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

SampleARRAY=($(for i in ${CandilesARRAY[@]}; do echo ${BARCODEID["${i}"]}; done | sort -u))


for i in ${SampleARRAY[@]}
  do

    echo -e " ********************** \n\n - Merging all ${i} samples and Re-Sorting - \n\n **********************"

    ls $OUT/${i}_*_sorted_rg.bam  > $OUT/${i}.bam.list
    samtools merge  -@ $THREADS -r $OUT/${i}_merged.bam -b $OUT/"${i}".bam.list
# removing the automatic deletion of intermediary bams to keep them in case something goes wrong
    # BAMARRAY=($(cat $OUT/${i}.bam.list))
    # for k in ${BAMARRAY[@]}
    #   do
    #     echo "Removing $k"
    #     rm $k
    # done

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
    -@ $THREADS -o $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam
    rm $OUT/${i}_cat_ref_sorted_rg_rmdup.bam
    samtools index $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam

    echo -e " ********************** \n\n - Realigning ${i} - \n\n **********************"

    # RealignerTargetCreator - ADDED -fixMisencodedQuals because these reads are from Illumina 1.5
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -nt $THREADS -R $REF -fixMisencodedQuals -I $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${i}_realignertargetcreator.intervals
    # IndelRealigner - ADDED -fixMisencodedQuals because these reads are from Illumina 1.5 and GATK was giving me an error.
    # Will try with whis and see if it works
    java -jar /home/tmp/Software/GATK_3.4/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $REF -fixMisencodedQuals -targetIntervals $OUT/${i}_realignertargetcreator.intervals \
    -I $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam \
    -o $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
    rm $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted.bam
    samtools flagstat $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    > $OUT/${i}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.stats

done
