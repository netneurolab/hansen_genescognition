# Molecular signatures of cognition and affect

This repository contains scripts, functions, and data I used or created in support of my work, ["Mapping gene transcription and neurocognition across human neocortex"](https://www.nature.com/articles/s41562-021-01082-z). Find the preprint [here](https://www.biorxiv.org/content/10.1101/2020.07.16.203026v1), the postprint [here](hansen2021postprint.pdf), and the article [here](hansen2021nathumbehav.pdf).
All analyses were run on Matlab version 9.8.0.1359463 (R2020a) Update 1.

Data and scripts are organized into five subfolders. The data that is used in multiple scripts is included in the root folder. This data includes:
- [gene_expression.mat](gene_expression.mat): node x gene matrix of normalized expression levels, created using [abagen](https://github.com/rmarkello/abagen).
- [label.mat](label.mat): a gene x 1 list of gene names, which correspond to the genes in gene_expression.mat
- [genes.mat](genes.mat): a struct of the indices that refer to stable genes in terms of three resolutions (34, 57, and 111 left-hemisphere nodes).
- [neurosynth.mat](neurosynth.mat): a node x term matrix of probabilistic measures that certain terms are pubslished alongside certain brain regions
- [nodes.mat](nodes.mat): a struct of the indices that refer to the left hemisphere brain regions
- [result.mat](result.mat): the original PLS result that is used in all other analyses
- [result34.mat](result34.mat): the PLS result computed on the 34-node parcellation, used in [scpt_brainspan.m](BrainSpan/scpt_brainspan.m).
- [spins.mat](spins.mat): a node x 10000 matrix of rotated left hemisphere brain regions used in spin tests
- [terms.mat](terms.mat): a struct of term names for Neurosynth terms and BrainMap terms

## The Main Analysis

The folder PLS contains the script [scpt_genes_cog_pls.m](PLS/scpt_genes_cog_pls.m) which performs partial least squares analysis on gene expression and functional activation matrices.
It also contains the script [scpt_cca.m](PLS/scpt_cca.m) which performs canonical correlation analysis on gene expression and functional activation matrices, and compares the results to the PLS results.
The significance of the latent variables is assessed against a permutation test that accounts for spatial autocorrelation.
The correlation of PLS-derived scores is cross-validated using the function [fcn_crossval_pls_brain_obvs.m](fcn_crossval_pls_brain_obvs.m) which assigns nodes on a distance-based method to account for spatial autocorrelation.
The terms that contribute most to the first latent variable are extrated.
Finally, PLS-derived scores are distributed among three network classifications: the intrinsic (resting-state) networks, the Von Economo cytoarchitectonic classes, and the Mesulam classes of laminar differentiation.

The main PLS code (`pls_analysis.m`) can be found [here](http://pls.rotman-baycrest.on.ca/source/), under "Latest PLS Applications".

The data in the folder is:
- [coords.mat](PLS/coords.mat): (x,y,z) coordinates of brain regions
- [rsn.mat](PLS/rsn.mat): resting-state network affiliation of brain regions
- [rsn_names.mat](PLS/rsn_names.mat): resting-state network names
- [ve_names.mat](PLS/ve_names.mat): Von Economo network names
- [voneconomo.mat](PLS/voneconomo.mat): Von Economo network affiliation of brain regions
- [mesulam_mapping.csv](PLS/mesulam_mapping.csv): Mesulam class affiliation of brain regions

## Gene Set Enrichment Analysis

The folder GO contains the script [scpt_GO.m](GO/scpt_GO.m) which performs gene set enrichment analysis based on two PLS-defined gene sets.
Analyses were adapted from [this repository](https://github.com/benfulcher/GeneSetEnrichmentAnalysis) which also provides two necessary files which can be found [here](https://figshare.com/s/71fe1d9b2386ec05f421). 

The data in the folder is:
- [gene_entrez_ids.csv](GO/gene_entrez_ids.csv): contains the entrezID corresponding to a list of genes. This file is a cleaner version of the original probes.csv from the Allen Institute. It was created using [abagen](https://github.com/rmarkello/abagen).

## Cell-Type Deconvolution

The folder CTD contains the script [scpt_ctd.m](CTD/scpt_ctd.m) which determines the ratio of genes that are preferentially expressed in seven different cell types.
Significance is assessed against a null model of random gene sets.
Cell type deconvolution comes from work discribed in [this paper](https://www.nature.com/articles/s41467-020-17051-5), and the data (alongside much more) can also be found at [Jakob Seidlitz's repo](https://github.com/jms290/PolySyn_MSNs)
The folder also contains the function [fcn_ctd.m](CTD/fcn_ctd.m) which was written after publication, but is easier to use so I've included it here.

The data in the folder is:
- [celltypes_PSP.csv](CTD/celltypes_PSP.csv) : A list of gene names preferentially expressed in each of seven cell types.

## Individual Differences in Behaviour

The folder HCP contains the script [scpt_hcp.m](HCP/scpt_hcp.m) which uses cortical thickness and T1w/T2w maps from the Human Connectome Project (S1200 release) to relate the PLS-derived gene score pattern to individual differences in behaviour.
Original data can be downloaded from [here](https://db.humanconnectome.org/data/projects/HCP_1200).
Note that the script is written for all 1096 subjects with full fMRI runs, but in reality only 417 unrelated subjects were used in analyses.
Due to privacy policies, their subject indices are not included.
This section got cut during the review process, so the results are only included in the preprint.

The data in the folder is: 
- [hcp_smyl_all_125.mat](HCP/hcp_smyl_all_125.mat): T1w/T2w ratios for all 1096 subjects with full fMRI runs parcellated into 219 cortical regions
- [hcp_thi_all_125.mat](HCP/hcp_thi_all_125.mat): cortical thickness for all 1096 subjects with full fMRI runs parcellated into 219 cortical regions

## Molecular Signature across Development

The folder BrainSpan contains the script [scpt_brainspan.m](BrainSpan/scpt_brainspan.m) which replicates results using gene expression estimates from the [BrainSpan](https://www.brainspan.org/static/download.html) database. 
The script also tracks the gene expression-functional activation signature across human development.
Many thanks to Jake Vogel for organizing the data for comparability with AHBA.

The data in the folder is:
- [gene_expression_AHBA_harmonized.csv](BrainSpan/gene_expression_AHBA_harmonized.csv): gene by sample matrix of extimated gene expression levels. Only includes genes that are also reported in the original AHBA dataset.
- [gene_metadata_AHBA_harmonized.csv](BrainSpan/gene_metadata_AHBA_harmonized.csv): includes information for each gene including gene name and entrez ID
- [samples_metadata.csv](BrainSpan/samples_metadata.csv): includes information for each sample, including age and brain region
- [mapping.mat](BrainSpan/mapping.mat) : an index-based map that links the 34 left hemisphere Desikan Killiany regions to the 16 cortical regions included in BrainSpan

