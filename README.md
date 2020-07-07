# Molecular signatures of cognition and affect

This repository contains scripts, functions, and data I used or created in support of my work-in-progress (soon to be preprint).
All analyses were run on Matlab version 9.8.0.1359463 (R2020a) Update 1.

## The Main Analysis

The folder PLS contains the script scpt_genes_cog_pls.m which performs partial least squares analysis on gene expression and functional activation matrices.
The significance of the latent variables is assessed against a permutation test that accounts for spatial autocorrelation.
The correlation of PLS-derived scores is cross-validated using the function fcn_crossval_pls_brain_obvs.m which assigns nodes on a distance-based method to account for spatial autocorrelation.
The terms that contribute most to the first latent variable are extrated.
Finally, PLS-derived scores are distributed among three network classifications: the intrinsic (resting-state) networks, the Von Economo cytoarchitectonic classes, and the Mesulam classes of laminar differentiation.

The data in the folder is:
- coords.mat : (x,y,z) coordinates of brain regions
- rsn.mat : resting-state network affiliation of brain regions
- rsn_names.mat : resting-state network names
- ve_names.mat : Von Economo network names
- voneconomo.mat : Von Economo network affiliation of brain regions
- mesulam_mapping.csv : Mesulam class affiliation of brain regions

## Gene Set Enrichment Analysis

The folder GO contains the script scpt_GO.m which performs gene set enrichment analysis based on two PLS-defined gene sets.
Analyses were adapted from [this repository][https://github.com/benfulcher/GeneSetEnrichmentAnalysis] which also provides two necessary files which can be found [here][https://figshare.com/s/71fe1d9b2386ec05f421]. 

The data in the folder is:
- gene_entrez_ids.csv : contains the entrezID corresponding to a list of genes

## Cell-Type Deconvolution

The folder CTD contains the script scpt_ctd.m which determines the ratio of genes that are preferentially expressed in seven different cell types.
Significance is assessed against a null model of random gene sets.
Cell type deconvolution comes from work discribed in [this paper][https://www.nature.com/articles/s41467-020-17051-5]. 

The data in the folder is:
- celltypes_PSP.csv : A list of gene names preferentially expressed in each of seven cell types.

## Individual Differences in Behaviour

The folder HCP contains the script scpt_hcp.m which uses cortical thickness and T1w/T2w maps from the Human Connectome Project (S1200 release) to relate the PLS-derived gene score pattern to individual differences in behaviour.
Original data can be downloaded from [here][https://db.humanconnectome.org/data/projects/HCP_1200]. 

The data in the folder is: 
- hcp_smyl_all_125.mat : T1w/T2w ratios for all 1096 subjects with full fMRI runs parcellated into 219 cortical regions
- hcp_thi_all_125.mat : cortical thickness for all 1096 subjects with full fMRI runs parcellated into 219 cortical regions

## Molecular Signature across Development

The folder BrainSpan contains the script scpt_brainspan.m which replicates results using gene expression estimates from the [BrainSpan][https://www.brainspan.org/static/download.html] database. 
The script also tracks the gene expression-functional activation signature across human development.

The data in the folder is:
- gene_expression_AHBA_harmonized.csv : gene by sample matrix of extimated gene expression levels. Only includes genes that are also reported in the original AHBA dataset.
- gene_metadata_AHBA_harmonized.csv : includes information for each gene including gene name and entrez ID
- samples_metadata.csv : includes information for each sample, including age and brain region
- mapping.mat : an index-based map that links the 34 left hemisphere Desikan Killiany regions to the 16 cortical regions included in BrainSpan

