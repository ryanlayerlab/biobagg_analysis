# biobagg_analysis

## Publication Figures
To generate all figures, run:
```
bash pubfigures/generate_all_figures.sh
```

### Figure 1: Workflow
<details>
  
![workflow](pub_figures/horizontal_workflow.001.jpeg)<br>

</details>

### Figure 2: Family Data
```
python plotting/figure2_related.py \
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
```
<details>
  
![figure2](pub_figures/figure2_family.png)<br>

</details>

### Figure 3.1: Ancestry Data
```
python plotting/figure3_ancestry.py \
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
```
<details>
  
![figure3_distribution](pub_figures/figure3_distribution.png)<br>
![figure3_topK](pub_figures/figure3_topk.png)<br>

</details>

#### Figure 3.2: 1KG subopulations
```
python plotting/figure_sup_subpops.py \
    --ancestry data/1kg_info/1kg_ancestry.tsv \
    --k 20 \
    --colors pub_figures/colors.txt \
    --genosis data/subpop_counts/genosis_counts.tsv \
    --dst data/subpop_counts/dst_counts.tsv \
    --pihat data/subpop_counts/pihat_counts.tsv \
    --kin data/subpop_counts/kinship_counts.tsv \
    --png pub_figures/
```
<details>

|                                                           |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![genosis-subpop-heatmap](pub_figures/genosis_counts.png) | ![dst-subpop-heatmap](pub_figures/dst_counts.png)         |
| ![pi-hat-subpop-heatmap](pub_figures/pihat_counts.png)    | ![kinship-subpop-heatmap](pub_figures/kinship_counts.png) |

</details>

### Figure 3.3: CCPM Ancestries
```
python plotting/plot_ccpm.py \
    --ancestry data/ccpm_data/R168_ccpm_f3_PCs_genetic_similarity_19Jun2024.txt
    --ancestry_dir data/ccpm_data/ccpm_results/
    --png pub_figures/ccpm_ancestry.png
```

<details>

![ccpm_ancestry](pub_figures/ccpm_ancestry.png)

</details>

### Figure X: FST
```
python plotting/plot_subpop_fst.py \
    --output pub_figures/subpop_fst.png \
    --input data/genosis.fst.txt data/kinship.fst.txt data/dst.fst.txt data/pihat.fst.txt \
    --labels GenoSiS King-robust DST Pi-Hat \
    --height 4 \
    --width 9 
    
python plottying/plot_ccpm_fst.py \    
    --ccpm_fst data/ccpm_data/plink_fst_my_fst.txt.fst.summary
    --out_dir pub_figures/
```
<details>

![fst](pub_figures/subpop_fst.png)
![ccpm_fst](pub_figures/ccpm_fst.png)


</details>

### Figure X: density
```
python plotting/plot_density.py \
    --ancestry_file data/1kg_info/1kg_ancestry.tsv" \
    --density_dir data/density_1kg/" \
    --distance_dir data/distance_1kg/" \
    --colors pub_figures/colors.txt" \
    --out pub_figures/ \
```
<details>

![hist](pub_figures/density_by_ancestry_chrm8.png)
![r2](pub_figures/r2_density_chrm8_seg76.png)

</details>


### Figure 4: Timing
```
tbd
```

### Figure 5: Cohort Quality
```
python plotting/figure5_quality.py \
    --ancestry data/1kg_info/1kg_ancestry.tsv \
    --k 20 \
    --colors pub_figures/colors.txt \
    --quality_dir data/quality_data/ \
    --png_hist pub_figures/figure5_
```
<details>
  
![figure5_cohort_quality_super_hist](pub_figures/figure5_super_histogram.png)<br>
![figure5_cohort_quality_sub_hist](pub_figures/figure5_sub_histogram.png)<br>
![figure5_cohort_quality_heatmap](pub_figures/figure5_heatmap.png)<br>


</details>



### SUPPLEMENTAL


--------------------------------------------
_!! **I have done some reorganizing of the data files, so paths might be wrong below** !!_


## 1KG population distributions
Plot 1KG population density plots
```
python plotting/plot_d1kg_distributions.py\
    --trios 1kg_trio_data/GenoSiS_ \
    --ancestry data/1kg_info/1kg_ancestry.tsv \
    --summary data/1KG_pop_hits.txt \
    --hits data/1kg_chr1-22_top_hits-100.txt \
    --ped data/1kg_info/1kg_trios.txt \
    --out 1kg_plots/
```

<details>
<summary>Example PNG:</summary>
  
![1kg_pop_k20](1kg_plots/1KG_pop_distributions.png)<br>

</details>


## GenoSiS v deCODE
```
python plotting/genosis_v_deCODE_smooth_hist.py \
    --POP data/decode_POP.txt \
    --IBD data/decode_IBD.txt \
    --png doc/genosis_v_decode.png \
    --height 3.5
```

<details>
<summary>Example PNG:</summary>
  
![genosis v decode](doc/genosis_v_decode.png)<br>

</details>



## 1KG trios
Label relationships in the 1KG trios data
```
python src/get_relations.py\
    --ped samples.trios\
    --ancestry data/1kg_info/1kg_ancestry.tsv \
    --pop AFR 
```
Plot GenoSiS scores for the 1KG trios data
```
python plotting/plot_1kg_trios.py \
    --ped data/1kg_info/1kg_trios.txt \
    --hits data/1kg_chr1-22_top_hits-100.txt \
    --dist data/1kg_info/1kg_trios \
    --pop AFR \
    --ancestry data/1kg_info/1kg_ancestry.tsv \
    --color darkorange 
```
<details>
<summary>Example PNG:</summary>

| Pop. | k=20                                                    | plink pi-hat                                         |
|------|---------------------------------------------------------|------------------------------------------------------|
AFR |  ![AFR_k=20](1kg_trio_plots/1KG_trios_AFR_20of20.png) | ![AFR_pihat](1kg_trio_plots/1KG_trios_plink_AFR.png) |
AMR |  ![AMR_k=20](1kg_trio_plots/1KG_trios_AMR_20of20.png) | ![AMR_pihat](1kg_trio_plots/1KG_trios_plink_AMR.png) |
EAS | ![EAS_k=20](1kg_trio_plots/1KG_trios_EAS_20of20.png) | ![EAS_pihat](1kg_trio_plots/1KG_trios_plink_EAS.png) |
EUR |  ![EUR_k=20](1kg_trio_plots/1KG_trios_EUR_20of20.png) | ![EUR_pihat](1kg_trio_plots/1KG_trios_plink_EUR.png) |
SAS |  ![SAS_k=20](1kg_trio_plots/1KG_trios_SAS_20of20.png) | ![SAS_pihat](1kg_trio_plots/1KG_trios_plink_SAS.png) |

</details>

## 1KG GenoSiS vs. Plink (kristen)
Compare 3 plink measurements to GenoSiS scores
```
python src/compare_genosis_plink.py \
    --genosis data/TOP_HITS_20.txt \
    --plinkgenome data/plink-genome.genome \
    --plinkkin data/plink2.kin0 \
    --top_k 20 \
    --ancestry data/1kg_info/1kg_ancestry.tsv 
    --pop AFR \
    --color darkorange 
```
| Pop. | k=20                                                          | k=100                                                           |
|------|---------------------------------------------------------------|-----------------------------------------------------------------|
AFR | ![AFR_k=20](1kg_plink_plots/genosis_plink_compare_AFR_20.png) | ![AFR_k=100](1kg_plink_plots/genosis_plink_compare_AFR_100.png) |
AMR | ![AMR_k=20](1kg_plink_plots/genosis_plink_compare_AMR_20.png) | ![AMR_k=100](1kg_plink_plots/genosis_plink_compare_AMR_100.png) |
EAS | ![EAS_k=20](1kg_plink_plots/genosis_plink_compare_EAS_20.png)    | ![EAS_k=100](1kg_plink_plots/genosis_plink_compare_EAS_100.png) |
EUR | ![EUR_k=20](1kg_plink_plots/genosis_plink_compare_EUR_20.png)    | ![EUR_k=100](1kg_plink_plots/genosis_plink_compare_EUR_100.png) |
SAS | ![SAS_k=20](1kg_plink_plots/genosis_plink_compare_SAS_20.png)    | ![SAS_k=100](1kg_plink_plots/genosis_plink_compare_SAS_100.png) |

## Targeted search
Get the SVS score distro (excluding self hits) to justify a 1.0 max 
cutoff
```
python plotting/plot_svs_score_hist.py -i "data/svs_results_100/chrm6.segment*" -o doc/svs_score_top_100.hist.png -N 1000000 -F 100


python plotting/plot_svs_score_hist.py -i "data/svs_results_20/chrm6.segment*" -o doc/svs_score_top_20.hist.png -N 1000000 -F 100
```

```
```
<details>
<summary>Example PNG:</summary>
  
![k=100](doc/svs_score_top_100.hist.png)<br>
![k=20](doc/svs_score_top_20.hist.png)<br>

</details>



## Look at iLASH v our calls
```
python src/ilash_pval.py \
    --ilash data/iLASH/chrm15.match \
    --segment_bed data/segment_boundary.bed.gz \
    --hapsis 'data/svs_results_chrm15-20/chrm15.segment*' \
> data/ilash_v_hapsis.txt

cat data/ilash_v_hapsis.txt \
| awk '{print $4,$3;}' \
| python src/scatter.py \
    -o doc/ilash_v_hapsis.chr15.png \
    --fig_x 6 --fig_y 3 \
    --x_label "GenoSiS popcount" \
    --y_label "iLASH IBD"
```

<details>
<summary>Example PNG:</summary>
  
![iLASH v Genosis](doc/ilash_v_hapsis.chr15.png)<br>

</details>




## Look at how comon / rare iLASH IBD calls are
```
cat data/iLASH/chrm15.match \
| cut -f 5-7 \
| sort \
| uniq -c \
| awk '{print $1;}' \
| python src/hist.py \
    -o doc/iLASH.chr15.ac.hist.png \
    -l \
    -b 50 \
    -x "IBD allele count" \
    -y "Freq." \
    -t "iLASH 1KG Chromosome 15" \
    --width 6 \
    --height 3

cat data/iLASH/chrm15.match \
| cut -f 5-7 \
| sort \
| uniq -c \
| awk '{print $4-$3,$1;}' \
| python src/scatter.py \
    -o doc/iLASH.chr15.len_x_ac.scatter.png \
    --log_x \
    --log_y \
    --x_label "IBD length" \
    --y_label "Allele count" \
    --fig_x 6 \
    --fig_y 3
```

<details>
<summary>Example PNG:</summary>
  
![iLASH chrm15 ac](doc/iLASH.chr15.ac.hist.png)<br>
![iLASH chrm15 len x ac](doc/iLASH.chr15.len_x_ac.scatter.png)<br>

</details>



## Compare the % of top k that match the query's super population

```
python plotting/percent_in_group.py \
    --pop data/population.txt \
    --knn data/1kg_chr1-22_top_hits.txt \
    --plink data/plink2.kin.top_20.txt \
    --out_file doc/hapsis_v_plink_in_pop.png

python plotting/percent_in_group.py \
    --pop data/population.txt \
    --knn data/1kg_chr1-22_top_100_hits.txt \
    --plink data/plink2.kin.top_100.txt \
    --out_file doc/hapsis_v_plink_in_pop.top_100.png
```

<details>
<summary>Example PNG:</summary>
  
![top 20](doc/hapsis_v_plink_in_pop.png)<br>
![top 20](doc/hapsis_v_plink_in_subpop.png)<br>
![top 100](doc/hapsis_v_plink_in_pop.top_100.png)<br>

</details>



## Generate violin plots for 1K subpopulation data

Arguments:
- population.txt: file with 1KG population information
- top_hits.txt: file with KNN hits (formatted with Ryan's aggregation)
- sub_population_violin: directory for png files

Example Run:
```
python plotting/evaluate_ancestry.py \
  --pop data/population.txt \
  --knn data/1kg_chr1-22_top_hits.txt \
  --png population_violin_plots/
```

<details>
<summary>Example PNGs:</summary>
  
![AFR: Esan](population_violin_plots/ESN.png)<br>
![AMR: Puerto Rican](population_violin_plots/PUR.png)<br>
![EAS: Japanese](population_violin_plots/JPT.png)<br>
![EUR: Finnish](population_violin_plots/FIN.png)<br>
![SAS: Bengali](population_violin_plots/BEB.png)<br>

</details>

## Generate violin plots for deCODE data

Arguments:
- decode_POP.txt: file with pop count scores for decode data 
- decode_IBD.txt: file with IBD scores for decode data 
- decode.png: violoin plot

Example Run:
```
python plotting/decode.py \
  --POP data/decode_POP.txt \
  --IBD data/decode_IBD.txt \
  --png doc/decode.png \
```

<details>
<summary>Example PNGs:</summary>

![deCODE: relationships, pop](doc/decode_rel_gess.png)<br>
![deCODE: relationships, pop + IBD](doc/decode_rel_both.png)<br>
![deCODE: meioses, pop](doc/decode_mei_gess.png)<br>
![deCODE: meioses, pop + IBD](doc/decode_mei_both.png)<br>
</details>

## Generate search time plot for 1K

Arguments:
- search.log: contains the timing in seconds for each segment query, which includes all samples 
- 3202: the number of samples, so we can get a per sample time

Example Run:
```
python plotting/search_time.py 
  --num_samples 3202 
  --in_file data/search.log \
  --out_file time_plots/1kg_search_time_histo.png \
  --height 4 \
  --width 8
per sample per segment median (ms) 1.6923797626483446
per sample per segment mean (ms) 1.7311439951686394
per sample per segment stdev (ms) 0.2026710271379037
per sample mean run time (s) 6.0624662710805755
```
<details>
<summary>Example PNG:</summary>

![](doc/1kg_search_time_histo.png)

</details>

## Look at the distribution of SVS scores

Arguments:
- chrm15-20.scores.gz: a comppress files where each line is a segment and all the top scores found for that segment

To get the scores file we condense scores to just one with:
```
for f in $( ls svs_results_chrm15-20/chrm*knn); do
    scores=$(cat $f | grep "_" | grep -v Query | cut -f2 | paste -sd " " -)
    echo -e "$f\t$scores"
done
```

Example Run:
```
python plotting/svs_scores.py \
    --in_file data/svs_results_chrm15-20.scores.gz \
    --out_file doc/svs_results_chrm15-20.scores.png \
    --height 8 \
    --width 5
15 [(15, 116, 3.836806668951722), (15, 28, 3.858078428857078), (15, 18, 11.544884247367353), (15, 39, 21.053507923862096), (15, 38, 47.9898720677839)]
16 [(16, 12, 3.1613097526350127), (16, 44, 3.5656954299731387), (16, 74, 4.133413055066931), (16, 2, 4.20666641220193), (16, 121, 9.45365259706147)]
17 [(17, 113, 2.926318128540011), (17, 2, 3.043968168123173), (17, 120, 3.1055936278300367), (17, 84, 3.1133913746537143), (17, 121, 4.326603849286288)]
18 [(18, 21, 3.0060410454243276), (18, 28, 3.060823315008067), (18, 30, 3.4353850232836374), (18, 9, 3.5373145689126986), (18, 102, 4.424481229262819)]
19 [(19, 1, 3.265355917200741), (19, 93, 3.4199648634766127), (19, 83, 3.6791844582441913), (19, 85, 5.445205323261027), (19, 2, 6.075723558402503)]
20 [(20, 2, 3.065332156191743), (20, 3, 5.15034327873235), (20, 89, 5.2502671985309455), (20, 94, 6.143665913587223), (20, 108, 6.660910728333406)]
```

<details>
<summary>Example PNG:</summary>

![](doc/svs_results_chrm15-20.scores.png)

</details>


## Look at populaiton structure just from top k

Arguments:
- 1kg_chr1-22_top_hits.txt: gives the top hist for each query 
- igsr-1000 genomes 30x on grch38.tsv: mapping from id to populations


```
python plotting/top_hits_umap.py \
  --in_file data/1kg_chr1-22_top_hits.txt \
  --label_file data/1kg_info/1kg_ancestry.tsv \
  --out_file doc/top_hit_umap.png
```

<details>
<summary>Example PNG:</summary>

![](doc/top_hit_umap.png)

</details>


```
python plotting/top_hits_pca.py \
  --in_file data/1kg_chr1-22_top_hits-100.txt \
  --label_file data/1kg_info/1kg_ancestry.tsv \
  --out_file doc/top_hit_pca.png
```


<details>
<summary>Example PNG:</summary>

![](doc/top_hit_pca.png)

</details>

## Compare the similarity of different groups using different metrics

Quick look at the distribution of scores from plink:
```
cat data/plink2.kin0 \
| tail -n +2 \
| cut -f 6 \
| python plotting/hist.py \
    -o doc/plink2.kin0.hist.png 
    --height 2 \
    --width 4 \
    -y Freq \
    -x "Plink2 Kinship" \
    -b 50 \
    -l
```

<details>
<summary>Example PNG:</summary>

![](doc/plink2.kin0.hist.png)

</details>

```
cat data/plink-genome.genome \
| tail -n +2 \
| awk '{print $10;}' \
| python plotting/hist.py \
    -o doc/plink-genome.genome.pi_hat.png \
    --height 2 \
    --width 4 \
    -y Freq \
    -x "Plink PI_HAT" \
    -b 50 \
    -l
```

<details>
<summary>Example PNG:</summary>

![](doc/plink-genome.genome.pi_hat.png)

</details>

Make pair files so the plotting is a little more uniform:
```
python src/sum_ilash.py \
    --data_dir data/iLASH/ \
    --out_file data/iLASH.pairs.txt

cat data/plink-genome.genome \
| tail -n +2 \
| awk '{OFS="\t"; print $2,$4,$10;}' \
> data/plink-genome.pairs.txt

cat data/plink2.kin0 \
| tail -n +2 \
| awk '{OFS="\t"; print $1,$2,$6;}' \
> data/plink2.kin.pairs.txt

cat data/iLASH.pairs.txt | cut -f 3 \
| python plotting/hist.py \
    -o doc/iLASH.pairs.png \
    --height 2 \
    --width 4 \
    -y Freq \
    -x "iLASH IBD sum" \
    -b 50 \
    -l
```

<details>
<summary>Example PNG:</summary>

![](doc/iLASH.pairs.png)

</details>


Look at the similarity score distributions for different groups using different scores:

<details>
<summary>Commands</summary>

```
python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits-100.txt \
    --pairs_file data/plink2.kin.pairs.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/distros_plink2-kin.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.2 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "Plink kinship"

python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits-100.txt \
    --pairs_file data/plink-genome.pairs.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/distros_plink-genome.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.45 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "Plink PI_HAT"

python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits-100.txt \
    --pairs_file data/iLASH.pairs.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/distros_iLASH.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.45 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "iLASH IBD"
```

</details>

<details>
<summary>Results:</summary>

![](doc/distros_plink2-kin.png)

![](doc/distros_plink-genome.png)

![](doc/distros_iLASH.png)

</details>

And as violin plots

```
python plotting/plot_distro_as_violins.py \
    --topk data/1kg_chr1-22_top_hits.txt \
    --pairs_file data/plink-genome.pairs.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/distros_plink-genome.violin.png \
    --height 2 \
    --width 8 \
    --y_label "Plink genome score"

python plotting/plot_distro_as_violins.py \
    --topk data/1kg_chr1-22_top_hits.txt \
    --pairs_file data/plink2.kin.pairs.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/distros_plink2-kin.violin.png \
    --height 2 \
    --width 8 \
    --y_label "Plink2 kinship"
```

<details>
<summary>Results:</summary>

![](doc/distros_plink-pihat.violin.png)

![](doc/distros_plink2-kin.violin.png)

![](doc/distros_plink-dist.violin.png)

</details>

## Compare rank of top k by different metrics

<details>
<summary>Commands</summary>

```
plot_ranks.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --pairs_file data/plink2.kin.pairs.txt \
    --out_file doc/ranks_plink2-kin.png \
    --height 2 \
    --width 8 \
    --y_label "Plink kinship rank"

~/src/biobagg_analysis$ python plotting/plot_ranks.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --pairs_file data/plink-genome.pairs.txt \
    --out_file doc/ranks_plink-genome.png \
    --height 2 \
    --width 8 \
    --y_label "Plink PI_HAT rank"
```

</details>

<details>
<summary>Results:</summary>

![](doc/ranks_plink2-kin.png)

![](doc/ranks_plink-genome.png)

</details>


## Look at the popcount scores distribution

```
python plotting/plot_topk_score_distro.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/topk_score_distro.png  \
    --height 2 \
    --width 8 \
    --bins 100 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt

python plotting/plot_topk_score_distro.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --label_file data/1kg_info/1kg_ancestry.tsv \
    --out_file doc/topk_score_distro.log.png  \
    --height 2 \
    --width 8 \
    --bins 100 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --log
```

<details>
<summary>Results:</summary>

![](doc/topk_score_distro.png)

![](doc/topk_score_distro.log.png)

</details>

## Look at top hits by population using plink

First coverst pairs files to top N files
```
python src/make_pairs_top_k.py \
    --pairs_file data/plink-genome.pairs.txt \
    --out_file data/plink-genome.top_20.txt \
    --N 20

python src/make_pairs_top_k.py \
    --pairs_file data/plink2.kin.pairs.txt \
    --out_file data/plink2.kin.top_20.txt \
    --N 20
```

Now make the violin plots
```
mkdir plink_genome_population_violin_plots/
python plotting/evaluate_ancestry.py \
    --pop data/population.txt \
    --knn data/plink-genome.top_20.txt \
    --png plink_genome_population_violin_plots/

mkdir plink2_kin_population_violin_plots/
python plotting/evaluate_ancestry.py \
    --pop data/population.txt \
    --knn data/plink2.kin.top_20.txt \
    --png plink2_kin_population_violin_plots/
```

<details>
<summary>Results:</summary>


| Pop. |  BIOBAG | Plink genome | Plink2 kinship|
|------|---------|---------------|---------------|
|  ACB  | ![](population_violin_plots/ACB.png) |![](plink_genome_population_violin_plots/ACB.png) |![](plink2_kin_population_violin_plots/ACB.png) |
|  ASW  | ![](population_violin_plots/ASW.png) |![](plink_genome_population_violin_plots/ASW.png) |![](plink2_kin_population_violin_plots/ASW.png) |
|  BEB  | ![](population_violin_plots/BEB.png) |![](plink_genome_population_violin_plots/BEB.png) |![](plink2_kin_population_violin_plots/BEB.png) |
|  CDX  | ![](population_violin_plots/CDX.png) |![](plink_genome_population_violin_plots/CDX.png) |![](plink2_kin_population_violin_plots/CDX.png) |
|  CEU  | ![](population_violin_plots/CEU.png) |![](plink_genome_population_violin_plots/CEU.png) |![](plink2_kin_population_violin_plots/CEU.png) |
|  CHB  | ![](population_violin_plots/CHB.png) |![](plink_genome_population_violin_plots/CHB.png) |![](plink2_kin_population_violin_plots/CHB.png) |
|  CHS  | ![](population_violin_plots/CHS.png) |![](plink_genome_population_violin_plots/CHS.png) |![](plink2_kin_population_violin_plots/CHS.png) |
|  CLM  | ![](population_violin_plots/CLM.png) |![](plink_genome_population_violin_plots/CLM.png) |![](plink2_kin_population_violin_plots/CLM.png) |
|  ESN  | ![](population_violin_plots/ESN.png) |![](plink_genome_population_violin_plots/ESN.png) |![](plink2_kin_population_violin_plots/ESN.png) |
|  FIN  | ![](population_violin_plots/FIN.png) |![](plink_genome_population_violin_plots/FIN.png) |![](plink2_kin_population_violin_plots/FIN.png) |
|  GBR  | ![](population_violin_plots/GBR.png) |![](plink_genome_population_violin_plots/GBR.png) |![](plink2_kin_population_violin_plots/GBR.png) |
|  GIH  | ![](population_violin_plots/GIH.png) |![](plink_genome_population_violin_plots/GIH.png) |![](plink2_kin_population_violin_plots/GIH.png) |
|  GWD  | ![](population_violin_plots/GWD.png) |![](plink_genome_population_violin_plots/GWD.png) |![](plink2_kin_population_violin_plots/GWD.png) |
|  IBS  | ![](population_violin_plots/IBS.png) |![](plink_genome_population_violin_plots/IBS.png) |![](plink2_kin_population_violin_plots/IBS.png) |
|  ITU  | ![](population_violin_plots/ITU.png) |![](plink_genome_population_violin_plots/ITU.png) |![](plink2_kin_population_violin_plots/ITU.png) |
|  JPT  | ![](population_violin_plots/JPT.png) |![](plink_genome_population_violin_plots/JPT.png) |![](plink2_kin_population_violin_plots/JPT.png) |
|  KHV  | ![](population_violin_plots/KHV.png) |![](plink_genome_population_violin_plots/KHV.png) |![](plink2_kin_population_violin_plots/KHV.png) |
|  LWK  | ![](population_violin_plots/LWK.png) |![](plink_genome_population_violin_plots/LWK.png) |![](plink2_kin_population_violin_plots/LWK.png) |
|  MSL  | ![](population_violin_plots/MSL.png) |![](plink_genome_population_violin_plots/MSL.png) |![](plink2_kin_population_violin_plots/MSL.png) |
|  MXL  | ![](population_violin_plots/MXL.png) |![](plink_genome_population_violin_plots/MXL.png) |![](plink2_kin_population_violin_plots/MXL.png) |
|  PEL  | ![](population_violin_plots/PEL.png) |![](plink_genome_population_violin_plots/PEL.png) |![](plink2_kin_population_violin_plots/PEL.png) |
|  PJL  | ![](population_violin_plots/PJL.png) |![](plink_genome_population_violin_plots/PJL.png) |![](plink2_kin_population_violin_plots/PJL.png) |
|  PUR  | ![](population_violin_plots/PUR.png) |![](plink_genome_population_violin_plots/PUR.png) |![](plink2_kin_population_violin_plots/PUR.png) |
|  STU  | ![](population_violin_plots/STU.png) |![](plink_genome_population_violin_plots/STU.png) |![](plink2_kin_population_violin_plots/STU.png) |
|  TSI  | ![](population_violin_plots/TSI.png) |![](plink_genome_population_violin_plots/TSI.png) |![](plink2_kin_population_violin_plots/TSI.png) |
|  YRI  | ![](population_violin_plots/YRI.png) |![](plink_genome_population_violin_plots/YRI.png) |![](plink2_kin_population_violin_plots/YRI.png) |

</details>

## Look at the top k for one all segemetns in one chromosome

There are 3 versions of this plot
1. All points are grey
2. Blue match the population and orange do not
3. The top k are colored

```
python plotting/query_result.py \
    --svs_results data/svs_results_chrm18.HG00096.bed.gz \
    --chrm 18 \
    --target HG00096 \
    --out_file doc/HG00096.png 

python plotting/query_result.py \
    --svs_results data/svs_results_chrm18.HG00096.bed.gz \
    --chrm 18 \
    --target HG00096 \
    --out_file doc/HG00096_labels.png \
    --label_file data/1kg_info/1kg_ancestry.tsv

python plotting/query_result.py \
    --svs_results data/svs_results_chrm18.HG00096.bed.gz \
    --chrm 18 \
    --target HG00096  \
    --out_file doc/HG00096_top_5.png \
    --k 5
```

<details>
<summary>Results:</summary>

![](doc/HG00096.png)

![](doc/HG00096_labels.png)

![](doc/HG00096_top_5.png)

</details>


