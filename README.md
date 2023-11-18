# biobagg_analysis

### Generate violin plots for 1K subpopulation data

Arguments:
- population.txt: file with 1KG population information
- top_hits.txt: file with KNN hits (formatted with Ryan's aggregation)
- sub_population_violin: directory for png files

Example Run:
```
python evaluate_ancestry.py \
  --pop population.txt \
  --knn top_hits.txt \
  --png sub_population_violin
```

Example PNGs:
![Esan](https://github.com/ryanlayerlab/biobagg_analysis/tree/main/population_violin_plots/ESN.png)<br>
![Puerto Rican](https://github.com/ryanlayerlab/biobagg_analysis/tree/main/population_violin_plots/PUR.png)<br>
![Japanese](https://github.com/ryanlayerlab/biobagg_analysis/tree/main/population_violin_plots/JPT.png)<br>
![Finnish](https://github.com/ryanlayerlab/biobagg_analysis/tree/main/population_violin_plots/FIN.png)<br>
![Bengali](https://github.com/ryanlayerlab/biobagg_analysis/tree/main/population_violin_plots/BEB.png)<br>
