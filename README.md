# Datasets and Results for QR-STAR paper

This repository includes the datasets and scripts used in the following papers:

- Tabatabaee, Y., Roch, S. and Warnow, T. (2023) Statistically consistent rooting of species trees under the multispecies coalescent model. In International Conference on Research in Computational Molecular Biology (RECOMB 2023) (pp. 41-57). https://doi.org/10.1101/2022.10.26.513897
- Tabatabaee, Y., Roch, S. and Warnow, T. (2023) QR-STAR: A Polynomial-Time Statistically Consistent Method for Rooting Species Trees Under the Coalescent. Journal of Computational Biology, 30(11), pp.1146-1181. https://www.liebertpub.com/doi/10.1089/cmb.2023.0185

For experiments in this study, we studied a collection of simulated datasets with incomplete lineage sorting (ILS). We used the 100- and 200-taxon published datasets from [Zhang et al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and [Mirarab and Warnow (2015)](https://academic.oup.com/bioinformatics/article/31/12/i44/215524) and the avian and mammalian simulated datasets from [Mirarab et al. (2014)](https://www.science.org/doi/full/10.1126/science.1250463).

All data can be accessed from [this](https://drive.google.com/drive/folders/1eGrH1ejoxBqVC6DVIhuU5cwJqjSGQZiU?usp=sharing) Google Drive link. This repository includes the scripts and results from the experiments in the paper.

## Simulated datasets
### 101-taxon simulations
This dataset, which was used for training and for designing QR-STAR, has four model conditions with varying sequence lengths (1600bp, 800bp, 400bp, 200bp) corresponding to different levels of gene tree estimation error (23%, 31%, 42%, and 55%), each with 50 replicates. The original dataset is from [Zhang et al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and available at [https://gitlab.com/esayyari/ASTRALIII/](https://gitlab.com/esayyari/ASTRALIII/). Results and intermediate data from the experiments in the paper are in `ASTRALIII.tar.xz`. Below is a description of files in each directory.

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
- `qr_le_v1.2.4_estgenetrees_[SEQ-LEN].non_[SHAPE-COEF]_[AB-RATIO]`: QR-STAR rooting of the true species tree given 1000 gene trees estimated from alignments of length `[SEQ-LEN]` where `[SHAPE-COEF]` and `[AB-RATIO]` specify the shape coefficient $C$ (default 1E-02) and the ratio $\frac{\alpha_{max}}{\beta_{min}}$ (default 0) in QR-STAR.

### 201-taxon simulations
This dataset, which was used for the testing experiments, has model conditions characterized by two different speciation rates and three tree heights (thus six tree shapes), and three number of genes for each tree shape. The number of replicates for each model condition is 50. The model conditions are named as `model.[NUM-SPECIES].[TREE-HEIGHT].[SPECIATION-RATE]` where `[NUM-SPECIES]` varies between 10, 50, 100, 200, 500 and 1000, `[SPECIATION-RATE]` varies between 1E-06 or 1E-07, and `[TREE-HEIGHT]` varies between 500K, 2M and 10M. The original dataset is from [Mirarab and Warnow (2015)](https://academic.oup.com/bioinformatics/article/31/12/i44/215524) and available at [https://datadryad.org/dataset/doi:10.6076/D10C7C](https://datadryad.org/dataset/doi:10.6076/D10C7C). Results and intermediate data from the experiments in the paper are in `ASTRALII.tar.xz`. Below is a description of files in each directory.

- `s_tree.trees`: true species tree in generation units
- `truegenetrees_[NUM-GENES]`: `[NUM-GENES]` true gene trees (by default 1000 if not specified)
- `estimatedgenetre_[NUM-GENES]`: `[NUM-GENES]` estimated gene trees (by default 1000 if not specified)
- `ad.txt`: average RF distance between the model species tree and true gene trees
- `gtee.txt`: average RF distance between estimated gene trees and true gene trees
- `astral.5.7.8-[GENE-TREES]_[NUM-GENES]`: ASTRAL species tree estimated using `[NUM-GENES]` estimated or true gene trees
- `rf_astral_[NUM-GENES].txt`: RF distance between the ASTRAL species tree estimated using `[NUM-GENES]` estimated gene trees and the true species tree
- `optimal_rooting_astral.5.7.8-[GENE-TREES]_[NUM-GENES]`: Optimal rooting of ASTRAL species tree estimated using `[NUM-GENES]` estimated or true gene trees
- `qr-v1.2.4-le-star.[GENE-TREES]_[NUM-GENES].[S-TREE]`: QR-STAR rooting of ASTRAL or true species tree (specified with `[S-TREE]`) using `[NUM-GENES]` estimated or true gene trees
- `qr-v1.2.4-le.[GENE-TREES]_[NUM-GENES].[S-TREE]`: QR rooting of ASTRAL or true species tree (specified with `[S-TREE]`) using `[NUM-GENES]` estimated or true gene trees

### Avian-like and mammalian-like simulations
These simulations use the 48-taxon avian-like and 37-taxon mammalian-like simulated datasets from [Mirarab et al. (2014)](https://www.science.org/doi/full/10.1126/science.1250463), which have model species trees based on biological datasets from [Jarvis et al. (2014)](https://www.science.org/doi/10.1126/science.1253451) and [Song et al. (2012)](https://www.pnas.org/doi/full/10.1073/pnas.1211733109), respectively. The default model condition (shown with `1X` ILS) has an ILS level that resembles the gene tree discordance in the corresponding biological data, but additional model conditions are created by multiplying or dividing branch lengths by two, thus decreasing or increasing the level of ILS, respectively (i.e., the highest ILS level is indicated by `0.5X`). The model conditions are named as `[ILS-level]-[NUM-GENES]-[SEQ-LEN]`, and there are 20 replicates in each condition. Original dataset is available from [https://doi.org/doi:10.5061/dryad.ht76hdrp0](https://doi.org/doi:10.5061/dryad.ht76hdrp0). Results and intermediate data from the experiments in the paper are in `avian-simulated.zip` and `mammalian-simulated.zip`. The model species trees are identical across all replicates and are named `avian-model-species.tre` and `mammalian-model-species.tre`. 

Below is a description of files in each directory.

- `gtee.txt`: average RF distance between estimated gene trees and true gene trees
- `estimatedgenetre_[NUM-GENES]`: `[NUM-GENES]` estimated gene trees (default value is 1000)
- `astral.5.7.8-estimatedgenetre_[NUM-GENES]`: ASTRAL species tree estimated using `[NUM-GENES]` estimated gene trees
- `optimal_rooting_astral.5.7.8-estimatedgenetre_[NUM-GENES]`: Optimal rooting of ASTRAL species tree estimated using `[NUM-GENES]` estimated gene trees
- `qr-v1.2.4-le-star.estimatedgenetre_[NUM-GENES].s_tree`: QR-STAR rooting of the true species tree given `[NUM-GENES]` estimated gene trees
- `qr-v1.2.4-le.estimatedgenetre_[NUM-GENES].s_tree`: QR rooting of the true species tree given `[NUM-GENES]` estimated gene trees
- `qr-v1.2.4-le-star.estimatedgenetre_[NUM-GENES].astral`: QR-STAR rooting of the ASTRAL species tree given `[NUM-GENES]` estimated gene trees
- `qr-v1.2.4-le.estimatedgenetre_[NUM-GENES].astral`: QR rooting of the ASTRAL species tree given `[NUM-GENES]` estimated gene trees
- `rf_astral_[NUM-GENES].txt`: RF distance between the ASTRAL species tree estimated using `[NUM-GENES]` estimated gene trees and the true species tree 
