#################
# Preparations: #
#################

# First we need to associate each sample to its fastq file, in order to calculate
# the total number of reads for each sample (TOTAL_SEQ).
# I will first copy the fastq Arrays, paths and BARCODESIDs of each project.

# List of all Lynx canadiens 3 (Murphy) sample codes in PN2017:
LCA_3ARRAY=($(ls /GRUPOS/grupolince/PN2017/LCA_3/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2 | uniq))
# List of all Lynx rufus 30 (Murphy) sample codes in PN2017:
LRU_30ARRAY=($(ls /GRUPOS/grupolince/PN2017/LRU_30/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2 | uniq))
# List of all Lynx rufus (Janecka) sample codes in PN2017/Bobcat1:
Bobcat1ARRAY=($(ls /GRUPOS/grupolince/PN2017/Bobcat1 | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3 | uniq))
# List of all Lynx pardinus sample codes in Lypa23:
CandilesARRAY=($(ls /home/ebazzicalupo/fastqs/Enrico_moveLypa/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | uniq))
# List of all Lynx lynx sample codes in MACROGEN project:
MacroGenARRAY=($(ls /backup/grupolince/raw_data/MACROGEN/MACROGEN_trimmed/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1 | uniq))
# List of all Lynx pardinus sample codes in Project Genoma part2:
PGenoma2ARRAY=($(ls /home/ebazzicalupo/fastqs/Enrico_move2/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3 | uniq))
# List of all Lynx pardinus sample codes in Project Genoma:
PGenomaARRAY=($(ls /home/ebazzicalupo/fastqs/Enrico_move/*.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3 | uniq))

# path to fastq files:
LCA_3PATH=/GRUPOS/grupolince/PN2017/LCA_3
LRU_30PATH=/GRUPOS/grupolince/PN2017/LRU_30
Bobcat1PATH=/GRUPOS/grupolince/PN2017/Bobcat1
CandilesPATH=/home/ebazzicalupo/fastqs/Enrico_moveLypa
MacroGenPATH=/backup/grupolince/raw_data/MACROGEN/MACROGEN_trimmed/
PGenoma2PATH=/home/ebazzicalupo/fastqs/Enrico_move2/
PGenomaPATH=/home/ebazzicalupo/fastqs/Enrico_move/

# BARCODES, where the ID of the fastq file is converted to our sample ID:
declare -A BARCODEID=(["LCA-3_1"]="c_lc_zz_0003" ["LCA-3_2"]="c_lc_zz_0003" ["LCA-3_3"]="c_lc_zz_0003" ["LCA-3_4"]="c_lc_zz_0003" ["LCA-3_5"]="c_lc_zz_0003" ["LCA-3_6"]="c_lc_zz_0003" ["LCA-3_7"]="c_lc_zz_0003" ["LCA-3_8"]="c_lc_zz_0003" ["LRU-30_1"]="c_lr_fl_0005" ["LRU-30_2"]="c_lr_fl_0005" ["LRU-30_3"]="c_lr_fl_0005" ["LRU-30_4"]="c_lr_fl_0005" ["LRU-30_5"]="c_lr_fl_0005" ["LRU-30_6"]="c_lr_fl_0005" ["LRU-30_7"]="c_lr_fl_0005" ["LRU-30_8"]="c_lr_fl_0005")
BARCODEID+=(["Bobcat1_S1_L001"]="c_lr_nm_0006" ["Bobcat1_S1_L002"]="c_lr_nm_0006" ["Bobcat1_S1_L003"]="c_lr_nm_0006" ["Bobcat1_S1_L004"]="c_lr_nm_0006")
BARCODEID+=(["6220RAAXX_lane3_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane4_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane6_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane7_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane8_sequence_0"]="c_lp_sm_0221" ["62AHEAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["621CYAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane5_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane2_sequence_0"]="c_lp_sm_0221")
BARCODEID+=(["LC1"]="c_lc_zz_0001" ["LL112"]="c_ll_vl_0112" ["LL146"]="c_ll_ya_0146" ["LL212"]="c_ll_cr_0212" ["LL90"]="c_ll_ki_0090" ["LR1"]="c_lr_zz_0001")
BARCODEID+=(["B09HCABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_5_0"]="c_lp_do_0173" ["B09HCABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_5_0"]="c_lp_do_0173" ["D0D6JABXX_4_0"]="c_lp_do_0443" ["D0D6JABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_4_0"]="c_lp_do_0443" ["B09HCABXX_3_0"]="c_lp_sm_0138" ["B09HCABXX_4_0"]="c_lp_sm_0138" ["B0B5KABXX_3_0"]="c_lp_sm_0138" ["B0B5KABXX_4_0"]="c_lp_sm_0138" ["C02CHABXX_1_0"]="c_lp_sm_0140" ["C02CHABXX_2_0"]="c_lp_sm_0140" ["C02CHABXX_4_0"]="c_lp_sm_0140" ["C02CHABXX_3_0"]="c_lp_sm_0140" ["D0D6JABXX_2_0"]="c_lp_sm_0185" ["D0D6JABXX_1_0"]="c_lp_sm_0185" ["B0999ABXX_2_0"]="c_lp_sm_0185" ["B0999ABXX_1_0"]="c_lp_sm_0185" ["D0D6JABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_8_0"]="c_lp_sm_0298" ["D0D6JABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_8_0"]="c_lp_sm_0298" ["C02CHABXX_6_0"]="c_lp_sm_0359" ["C02CHABXX_7_0"]="c_lp_sm_0359" ["C02CHABXX_5_0"]="c_lp_sm_0359" ["C02CHABXX_8_0"]="c_lp_sm_0359" ["B09HCABXX_8_0"]="h_lp_do_0007" ["B09HCABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_8_0"]="h_lp_do_0007")

# Then a small loop for each project will create a file for each fastqs that contains
# the count of the number of lines of the fastq devided by 4 (=number of reads in that fastq).
# These files will be used to later calculate the sum of all the number of reads from all of the
# fastqs of each sample (column 2 of the table = total_seq).

for i in ${LCA_3ARRAY[@]}
  do
    echo $i
    zcat $LCA_3PATH/${i}_AGTTCC_R1.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $LCA_3PATH/${i}_AGTTCC_R2.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${LRU_30ARRAY[@]}
  do
    echo $i
    zcat $LRU_30PATH/${i}_GGTAGC_R1.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $LRU_30PATH/${i}_GGTAGC_R2.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${Bobcat1ARRAY[@]}
  do
    echo $i
    zcat $Bobcat1PATH/${i}_R1_001.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $Bobcat1PATH/${i}_R2_001.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${CandilesARRAY[@]}
  do
    echo $i
    zcat $CandilesPATH/${i}_1.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $CandilesPATH/${i}_2.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${MacroGenARRAY[@]}
  do
    echo $i
    zcat $MacroGenPATH/${i}_R1_trimmed.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $MacroGenPATH/${i}_R2_trimmed.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${PGenoma2ARRAY[@]}
  do
    echo $i
    zcat $PGenoma2PATH/${i}_1.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $PGenoma2PATH/${i}_2.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done

for i in ${PGenomaARRAY[@]}
  do
    echo $i
    zcat $PGenomaPATH/${i}_1.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar1.rawseq
    zcat $PGenomaPATH/${i}_2.fastq.gz | wc -l | awk '{print $1/4}' > /home/ebazzicalupo/mapstats/${BARCODEID["${i}"]}_${i}.borrar2.rawseq
done


############################
# Now to the actual table! #
############################

# List of all the BAMS
SAMPLELIST=($(ls ~/CatRef_bams/*_indelrealigner.bam | rev | cut -d'/' -f 1 | rev | cut -d "_" -f1-4 | sort | uniq ))

# Create file with table header
echo "sample_name,total_seq,total_reads,duplicates,mapped,properly_pair,\
with_mat_to_another_chr,uniq_mapped,uniq_mapped_vs_total_reads,duplicates_vs_total_reads,\
mapped_vs_total_seq,properly_pair_vs_total_seq,with_mat_to_another_chr_vs_total_seq,\
samtools_depth,stdev_samtools_depth,cov_gt_zero,cov_gt_five,cov_gt_twenty,cov_gt_fifty,cov_gt_seventy" \
> ~/mapstats/mapstats.csv

# Loop adding a row for each sample
for sample in "${SAMPLELIST[@]}"
  do

    echo "${sample}"

    # column 1 : sample_name --> self explanatory
    NAME="${sample}"

    # column 2 : total_seq --> it's the number of reads in the fastq file
    TOTAL_SEQ="$(cat /home/ebazzicalupo/mapstats/"${sample}"*.rawseq | awk '{sum+=$1}END{print sum}')"

    # column 3 : total_reads --> number of alignments in the BAM file
    TOTAL_READS="$(grep "in total" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 4 : duplicates --> marked duplicates
    DUPLICATES="$(grep "duplicates" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 5 : mapped --> mapped reads
    MAPPED="$(grep "0 mapped" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 6 : properly_pair --> pairs of reads properly paired
    PROPERLY_PAIR="$(grep "properly paired" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 7 : with_mat_to_another_chr --> reads with mate mapped to another chromosome
    WITH_MATE_TO_ANOTHER_CHR="$(grep "with mate mapped to a different chr$" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 8 : uniq_mapped --> mapped reads without duplicates
    UNIQ_MAPPED="$(echo "$MAPPED - $DUPLICATES" | bc )"

    # column 9 : uniq_mapped_vs_total_reads
    UNIQ_MAPPED_vs_TOTAL_READS="$(echo "scale=5; $UNIQ_MAPPED *100 / $TOTAL_READS"  | bc -l)"

    # column 10 : duplicates_vs_total_reads
    DUPLICATES_vs_TOTAL_READS="$(echo "scale=5; $DUPLICATES *100 / $TOTAL_READS" | bc -l )"

    # column 11 : MAPPED_vs_TOTAL_SEQ
    MAPPED_vs_TOTAL_SEQ="$(echo "scale=5; $MAPPED *100  / $TOTAL_SEQ" | bc )"

    # column 12 : properly_pair_vs_total_seq
    PROPERLY_PAIR_vs_TOTAL_SEQ="$(echo "scale=5; $PROPERLY_PAIR *100 / $TOTAL_SEQ" | bc )"

    # column 13 : with_mat_to_another_chr_vs_total_seq
    WITH_MATE_TO_ANOTHER_CHR_vs_TOTAL_SEQ="$(echo "scale=5; $WITH_MATE_TO_ANOTHER_CHR *100 / $TOTAL_SEQ" | bc )"

    # column 14 : samtools_depth --> mean depth calculated with samtools depth (whole reference genome length = 2563897203)
    SAMTOOLS_DEPTH="$(samtools depth ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/2563897203}')"

    # column 15 : stdev_samtools_depth --> standard deviation of mean depth
    STDEV_SAMTOOLS_DEPTH="$(samtools depth ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sqrt(sumsq/2563897203 - (sum/2563897203)**2)}')"

    # column 16 : cov_gt_zero --> percentage of ref. genome covered by at least 1 read
    COV_GT_ZERO="$(samtools mpileup ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=1 '$4>=X' | wc -l)"
    COV_GT_ZERO_PC="$(echo "scale=5; $COV_GT_ZERO *100 / 2563897203" | bc )"

    # column 17 : cov_gt_five --> percentage of ref. genome covered by at least 5 read
    COV_GT_FIVE="$(samtools mpileup ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=5 '$4>=X' | wc -l)"
    COV_GT_FIVE_PC="$(echo "scale=5; $COV_GT_FIVE *100 / 2563897203" | bc )"

    # column 18 : cov_gt_twenty --> percentage of ref. genome covered by at least 20 read
    COV_GT_TWENTY="$(samtools mpileup ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=20 '$4>=X' | wc -l)"
    COV_GT_TWENTY_PC="$(echo "scale=5; $COV_GT_TWENTY *100 / 2563897203" | bc )"

    # column 19 : cov_gt_fifty --> percentage of ref. genome covered by at least 50 read
    COV_GT_FIFTY="$(samtools mpileup ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=50 '$4>=X' | wc -l)"
    COV_GT_FIFTY_PC="$(echo "scale=5; $COV_GT_FIFTY *100 / 2563897203" | bc )"

    # column 20 : cov_gt_fifty --> percentage of ref. genome covered by at least 70 read
    COV_GT_SEVENTY="$(samtools mpileup ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=70 '$4>=X' | wc -l)"
    COV_GT_SEVENTY_PC="$(echo "scale=5; $COV_GT_SEVENTY *100 / 2563897203" | bc )"

    echo "$NAME,$TOTAL_SEQ,$TOTAL_READS,$DUPLICATES,$MAPPED,$PROPERLY_PAIR,$WITH_MATE_TO_ANOTHER_CHR,$UNIQ_MAPPED,$UNIQ_MAPPED_vs_TOTAL_READS,$DUPLICATES_vs_TOTAL_READS,$MAPPED_vs_TOTAL_SEQ,$PROPERLY_PAIR_vs_TOTAL_SEQ,$WITH_MATE_TO_ANOTHER_CHR_vs_TOTAL_SEQ,$SAMTOOLS_DEPTH,$STDEV_SAMTOOLS_DEPTH,$COV_GT_ZERO_PC,$COV_GT_FIVE_PC,$COV_GT_TWENTY_PC,$COV_GT_FIFTY_PC,$COV_GT_SEVENTY_PC" >> ~/mapstats/mapstats.csv

done

rm /home/ebazzicalupo/mapstats/*borrar*
