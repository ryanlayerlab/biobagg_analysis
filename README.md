# biobagg_analysis

## Generate violin plots for 1K subpopulation data

Arguments:
- population.txt: file with 1KG population information
- top_hits.txt: file with KNN hits (formatted with Ryan's aggregation)
- sub_population_violin: directory for png files

Example Run:
```
python plotting/evaluate_ancestry.py \
  --pop population.txt \
  --knn top_hits.txt \
  --png population_violin_plots
```

<details>
<summary>Example PNGs:</summary>
  
![AFR: Esan](population_violin_plots/ESN.png)<br>
![AMR: Puerto Rican](population_violin_plots/PUR.png)<br>
![EAS: Japanese](population_violin_plots/JPT.png)<br>
![EUR: Finnish](population_violin_plots/FIN.png)<br>
![SAS: Bengali](population_violin_plots/BEB.png)<br>

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
  --label_file data/igsr-1000\ genomes\ 30x\ on\ grch38.tsv \
  --out_file doc/top_hit_umap.png
```

<details>
<summary>Example PNG:</summary>

![](doc/top_hit_umap.png)

</details>


```
python plotting/top_hits_pca.py \
  --in_file data/1kg_chr1-22_top_hits.txt \
  --label_file data/igsr-1000\ genomes\ 30x\ on\ grch38.tsv \
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

```
python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --pairs_file data/plink2.kin.pairs.txt \
    --label_file data/igsr-1000\ genomes\ 30x\ on\ grch38.tsv \
    --out_file doc/distros_plink2-kin.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.2 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "Plink kinship"

python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --pairs_file data/plink-genome.pairs.txt \
    --label_file data/igsr-1000\ genomes\ 30x\ on\ grch38.tsv \
    --out_file doc/distros_plink-genome.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.45 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "Plink PI_HAT"

python plotting/plot_distros.py \
    --topk_file data/1kg_chr1-22_top_hits.txt \
    --pairs_file data/iLASH.pairs.txt \
    --label_file data/igsr-1000\ genomes\ 30x\ on\ grch38.tsv \
    --out_file doc/distros_iLASH.png \
    --width 2 \
    --height 3 \
    --bins 50 \
    --alpha 1 \
    --outlier 0.45 \
    --ped_file data/1kGP.3202_samples.pedigree_info.txt \
    --x_label "iLASH IBD"
```
<details>
<summary>Example PNG:</summary>

![](doc/distros_plink2-kin.png)

![](doc/distros_plink-genome.png)

![](doc/distros_iLASH.png)

</details>




