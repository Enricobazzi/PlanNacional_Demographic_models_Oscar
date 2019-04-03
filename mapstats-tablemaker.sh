
# List of all the samples
SAMPLELIST=($(ls ~/CatRef_bams/*_indelrealigner.bam | rev | cut -d'/' -f 1 | rev | cut -d "_" -f1-4 | sort | uniq ))

# Create file with table header
echo "sample_name,total_seq,total_reads,duplicates,mapped,properly_pair,\
with_mat_to_another_chr,uniq_mapped,uniq_mapped_vs_total_reads,duplicates_vs_total_reads,\
mapped_vs_total_seq,properly_pair_vs_total_seq,with_mat_to_another_chr_vs_total_seq,\
samtools_depth,stdev_samtools_depth,cov_gt_zero,cov_gt_five,cov_gt_twenty,cov_gt_fifty" \
> raw.stats.csv

# Loop adding a row for each sample
for sample in "${SAMPLELIST[@]}"
  do

    echo "${sample}"

    # column 1 : sample_name
    NAME="${sample}"

    # column 2 : total_seq - Still have to figure out how this is calculated!
    TOTAL_SEQ="$(cat /home/GRUPOS/grupolince/lynx_genomes_5x/BAM_files_final/stats_information/"${sample}"*.rawseq | awk '{sum+=$1}END{print sum}')"

    # column 3 : total_reads
    TOTAL_READS="$(grep "in total" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 4 : duplicates
    DUPLICATES="$(grep "duplicates" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 5 : mapped
    MAPPED="$(grep "0 mapped" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 6 : properly_pair
    PROPERLY_PAIR="$(grep "properly paired" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 7 : with_mat_to_another_chr
    WITH_MATE_TO_ANOTHER_CHR="$(grep "with mate mapped to a different chr$" ~/CatRef_bams/"${sample}"*.stats | cut -d "+" -f 1)"

    # column 8 : uniq_mapped
    UNIQ_MAPPED="$(echo "$MAPPED - $DUPLICATES" | bc )"

    # column 9 : uniq_mapped_vs_total_reads
    UNIQ_MAPPED_vs_TOTAL_READS="$(echo "scale=5; $UNIQ_MAPPED *100 / $TOTAL_READS"  | bc -l)"

    # column 10 : duplicates_vs_total_reads
    DUPLICATES_vs_TOTAL_READS="$(echo "scale=5; $DUPLICATES / $TOTAL_READS" | bc -l )"

    # column 11 : MAPPED_vs_TOTAL_SEQ
    MAPPED_vs_TOTAL_SEQ="$(echo "scale=5; $MAPPED *100  / $TOTAL_SEQ" | bc )"

    # column 12 : properly_pair_vs_total_seq
    PROPERLY_PAIR_vs_TOTAL_SEQ="$(echo "scale=5; $PROPERLY_PAIR *100 / $TOTAL_SEQ" | bc )"

    # column 13 : with_mat_to_another_chr_vs_total_seq
    WITH_MATE_TO_ANOTHER_CHR_vs_TOTAL_SEQ="$(echo "scale=5; $WITH_MATE_TO_ANOTHER_CHR *100 / $TOTAL_SEQ" | bc )"

    # column 14 : samtools_depth
    SAMTOOLS_DEPTH="$(samtools depth ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/2563897203}')"

    # column 15 : stdev_samtools_depth
    STDEV_SAMTOOLS_DEPTH="$(samtools depth ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sqrt(sumsq/2563897203 - (sum/2563897203)**2)}')"

    # column 16 : cov_gt_zero
    COV_GT_ZERO="$(samtools ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=1 '$4>=X' | wc -l)"

    # column 17 : cov_gt_five
    COV_GT_FIVE="$(samtools ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=5 '$4>=X' | wc -l)"

    # column 18 : cov_gt_twenty
    COV_GT_TWENTY="$(samtools ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=20 '$4>=X' | wc -l)"

    # column 19 : cov_gt_fifty
    COV_GT_FIFTY="$(samtools ~/CatRef_bams/"${sample}"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | awk -v X=50 '$4>=X' | wc -l)"

done
