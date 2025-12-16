# Data and Results for QR-STAR

This repository includes the datasets and scripts used in the following papers:

- Tabatabaee, Y., Roch, S. and Warnow, T. (2023) Statistically consistent rooting of species trees under the multispecies coalescent model. In International Conference on Research in Computational Molecular Biology (RECOMB 2023) (pp. 41-57). https://doi.org/10.1101/2022.10.26.513897
- Tabatabaee, Y., Roch, S. and Warnow, T. (2023) QR-STAR: A Polynomial-Time Statistically Consistent Method for Rooting Species Trees Under the Coalescent. Journal of Computational Biology, 30(11), pp.1146-1181. https://www.liebertpub.com/doi/10.1089/cmb.2023.0185

For experiments in this study, we studied a collection of simulated datasets with incomplete lineage sorting (ILS). We used the 100- and 200-taxon published datasets from [Zhang et. al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and [Mirarab and Warnow (2015)](https://academic.oup.com/bioinformatics/article/31/12/i44/215524) and the avian and mammalian simulated datasets from [Mirarab et al. (2014)](https://www.science.org/doi/full/10.1126/science.1250463).

All data can be accessed from [this](https://drive.google.com/drive/folders/1eGrH1ejoxBqVC6DVIhuU5cwJqjSGQZiU?usp=sharing) Google Drive link. This repository includes the scripts and results from the experiments in the paper.

## Simulated datasets
### 100-taxon simulations
This dataset, which was used for training and for designing QR-STAR, has four model conditions with varying sequence lengths (1600bp, 800bp, 400bp, 200bp) corresponding to different levels of gene tree estimation error (23%, 31%, 42%, and 55%), each with 50 replicates. The original dataset is from [Zhang et. al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and available at [https://gitlab.com/esayyari/ASTRALIII/](https://gitlab.com/esayyari/ASTRALIII/). Results and intermediate data from the experiments in the paper are in `ASTRALIII.tar.xz`. Below is a description of files in each directory.

- `s_tree.trees`: true species tree in generation units
- `truegenetrees_[NUM-GENES]`: `[NUM-GENES]` true gene trees (by default 1000 if not specified)
- `estimatedgenetre_[SEQ-LEN].gtr.rerooted.final.contracted.non_[NUM-GENES]`: `[NUM-GENES]` gene trees estimated from alignments of length `[SEQ-LEN]`
- `ad.txt`: average RF distance between the model species tree and true gene trees
- `gtee_[SEQ-LEN].txt`: average RF distance between estimated gene trees (inferred from alignments of length `[SEQ-LEN]`) and true gene trees
- `astral.5.7.8-estimatedgenetre_[SEQ-LEN].non_[NUM-GENES]`: ASTRAL species tree estimated using `[NUM-GENES]` gene trees estimated from alignments of length `[SEQ-LEN]`
- `astral.5.7.8-truegenetrees_[NUM-GENES]`: ASTRAL species tree estimated using `[NUM-GENES]` true gene trees
- `optimal_rooting_astral.5.7.8-estimatedgenetre_[SEQ-LEN].non_[NUM-GENES]`: Optimal rooting of ASTRAL species tree estimated using `[NUM-GENES]` gene trees estimated from alignments of length `[SEQ-LEN]`
- `optimal_rooting_astral.5.7.8-truegenetrees_[NUM-GENES]`: Optimal rooting of ASTRAL species tree estimated using `[NUM-GENES]` true gene trees
- `qr-v1.2.4-le-star.estimatedgenetre_[SEQ-LEN].non_[NUM-GENES].astral`: QR-STAR rooting of ASTRAL species tree estimated using `[NUM-GENES]` gene trees estimated from alignments of length `[SEQ-LEN]`
- `qr-v1.2.4-le.estimatedgenetre_[SEQ-LEN].non_[NUM-GENES].astral`: QR rooting of ASTRAL species tree estimated using `[NUM-GENES]` gene trees estimated from alignments of length `[SEQ-LEN]`
- ``

### 200-taxon simulations

### Avian-like simulations

### Mammalian-like simulations
