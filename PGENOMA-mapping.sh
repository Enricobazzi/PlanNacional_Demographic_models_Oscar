#!/bin/bash

#####################################################################################
## Mapping GENOMA PROJECT FASTQ sequences (Illumina 1.9) to FELIS CATUS reference genome ##
#####################################################################################

# After the MACROGEN mapping, this script will be used to go through the same steps
# (mapping, adding read groups and merging, marking duplicates and realigning) with
# the samples of the GENOMA project.

# The MAIN DIFFERENCE with the Macrogen mapping will be that these samples have been
# sequenced in more than one lane, so the loop will be split in two parts.
# The first loop will generate a bam file for each run. The second lopp will
# merge the bams belonging to the same individual from the different runs and
# run the rest of the steps on the merged bam.

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

# List of all Lynx pardinus sample codes in Project Genoma:
PGenomaARRAY=($(ls /home/ebazzicalupo/fastqs/Enrico_move/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3 | uniq))
# Path to cat reference genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
# No. of computer cores used. 20 = OK, >20 = ask people first!
THREADS=10
# Path to output files, were BAMS are generated:
OUT=/home/ebazzicalupo/CatRef_bams
# path to Project Genoma fastq files:
PGenomaPATH=/home/ebazzicalupo/fastqs/Enrico_move/
# BARCODES PGenoma, where the ID given by the sequencing company is converted to our sample ID:
declare -A BARCODEID=(["B09HCABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_5_0"]="c_lp_do_0173" ["B09HCABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_5_0"]="c_lp_do_0173" ["D0D6JABXX_4_0"]="c_lp_do_0443" ["D0D6JABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_4_0"]="c_lp_do_0443" ["B09HCABXX_3_0"]="c_lp_sm_0138" ["B09HCABXX_4_0"]="c_lp_sm_0138" ["B0B5KABXX_3_0"]="c_lp_sm_0138" ["B0B5KABXX_4_0"]="c_lp_sm_0138" ["C02CHABXX_1_0"]="c_lp_sm_0140" ["C02CHABXX_2_0"]="c_lp_sm_0140" ["C02CHABXX_4_0"]="c_lp_sm_0140" ["C02CHABXX_3_0"]="c_lp_sm_0140" ["D0D6JABXX_2_0"]="c_lp_sm_0185" ["D0D6JABXX_1_0"]="c_lp_sm_0185" ["B0999ABXX_2_0"]="c_lp_sm_0185" ["B0999ABXX_1_0"]="c_lp_sm_0185" ["D0D6JABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_8_0"]="c_lp_sm_0298" ["D0D6JABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_8_0"]="c_lp_sm_0298" ["C02CHABXX_6_0"]="c_lp_sm_0359" ["C02CHABXX_7_0"]="c_lp_sm_0359" ["C02CHABXX_5_0"]="c_lp_sm_0359" ["C02CHABXX_8_0"]="c_lp_sm_0359" ["B09HCABXX_8_0"]="h_lp_do_0007" ["B09HCABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_8_0"]="h_lp_do_0007")


##########################
## Mapping with BWA MEM ##
##########################

for i in ${PGenomaARRAY[@]}
  do

    echo -e " ********************** \n\n - Mapping ${i} - \n\n **********************"
    bwa mem $REF $PGenomaPATH/${i}_1.fastq.gz $PGenomaPATH/${i}_2.fastq.gz \
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

SampleARRAY=($(ls $OUT/*_lp_*.bam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | sort -u))

for i in ${SampleARRAY[@]}
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
