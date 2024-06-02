# if any of the figures fail to generate, the script will stop
set -e pipefail

# figure 1 = workflow

# figure 2 = family data
echo 'Plotting figure 2...'
echo '...figure2_family.png'

python -W ignore plotting/figure2_related.py \
  --colors pub_figures/colors.txt \
  --decode_genosis data/decode/decode_POP.txt \
  --decode_ibd data/decode/decode_IBD.txt \
  --AFR_genosis data/1kg_trio_data/GenoSiS_AFR_trio_scores_20.txt \
  --AMR_genosis data/1kg_trio_data/GenoSiS_AMR_trio_scores_20.txt \
  --EAS_genosis data/1kg_trio_data/GenoSiS_EAS_trio_scores_20.txt \
  --EUR_genosis data/1kg_trio_data/GenoSiS_EUR_trio_scores_20.txt \
  --SAS_genosis data/1kg_trio_data/GenoSiS_SAS_trio_scores_20.txt \
  --AFR_dst data/1kg_trio_data/plink_dst_AFR_trio_scores_20.txt \
  --AMR_dst data/1kg_trio_data/plink_dst_AMR_trio_scores_20.txt \
  --EAS_dst data/1kg_trio_data/plink_dst_EAS_trio_scores_20.txt \
  --EUR_dst data/1kg_trio_data/plink_dst_EUR_trio_scores_20.txt \
  --SAS_dst data/1kg_trio_data/plink_dst_SAS_trio_scores_20.txt \
  --AFR_pihat data/1kg_trio_data/plink_pihat_AFR_trio_scores_20.txt \
  --AMR_pihat data/1kg_trio_data/plink_pihat_AMR_trio_scores_20.txt \
  --EAS_pihat data/1kg_trio_data/plink_pihat_EAS_trio_scores_20.txt \
  --EUR_pihat data/1kg_trio_data/plink_pihat_EUR_trio_scores_20.txt \
  --SAS_pihat data/1kg_trio_data/plink_pihat_SAS_trio_scores_20.txt \
  --AFR_kin data/1kg_trio_data/plink_kin_AFR_trio_scores_20.txt \
  --AMR_kin data/1kg_trio_data/plink_kin_AMR_trio_scores_20.txt \
  --EAS_kin data/1kg_trio_data/plink_kin_EAS_trio_scores_20.txt \
  --EUR_kin data/1kg_trio_data/plink_kin_EUR_trio_scores_20.txt \
  --SAS_kin data/1kg_trio_data/plink_kin_SAS_trio_scores_20.txt \
  --png pub_figures/figure2_family.png

# figure 3 = population data
echo 'Plotting figure 3...'
echo '...figure3_distribution.png'
echo '...figure3_topk.png'

python -W ignore plotting/figure3_ancestry.py \
  --ancestry data/1kg_info/1kg_ancestry.tsv \
  --k 20 \
  --colors pub_figures/colors.txt \
  --genosis_groups data/1kg_pop_hits.txt \
  --genosis_k data/1kg_top_hits/TOP_HITS_20.txt \
  --dst_groups data/1kg_plink_topK/plink_DST_20_groups.txt \
  --pihat_groups data/1kg_plink_topK/plink_pihat_20_groups.txt \
  --kinship_groups data/1kg_plink_topK/plink_kin_20_groups.txt \
  --dst_k data/1kg_plink_topK/plink_DST_top_20.txt \
  --pihat_k data/1kg_plink_topK/plink_pihat_top_20.txt \
  --kinship_k data/1kg_plink_topK/plink_kin_top_20.txt \
  --png_dist pub_figures/figure3_distribution.png \
  --png_k pub_figures/figure3_topk.png

# figure 4 = timing

# figure 5 = cohort quality
echo 'Plotting figure 5...'
echo '...figure5_quality.png'

python -W ignore plotting/figure5_quality.py \
  --ancestry data/1kg_info/1kg_ancestry.tsv \
  --k 20 \
  --colors pub_figures/colors.txt \
  --quality_dir data/quality_data/ \
  --png pub_figures/figure5_

# supplementary figures
echo 'Plotting supplementary figures...'
echo '...supplement_subpop_genosis.png'
echo '...supplement_subpop_dst.png'
echo '...supplement_subpop_pihat.png'
echo '...supplement_subpop_kinship.png'

python -W ignore plotting/figure_sup_subpops.py \
  --ancestry data/1kg_info/1kg_ancestry.tsv \
  --k 20 \
  --colors pub_figures/colors.txt \
  --genosis data/subpop_counts/genosis_subpop_counts.tsv \
  --dst data/subpop_counts/dst_subpop_counts.tsv \
  --pihat data/subpop_counts/pihat_subpop_counts.tsv \
  --kin data/subpop_counts/kinship_subpop_counts.tsv \
  --png pub_figures/

