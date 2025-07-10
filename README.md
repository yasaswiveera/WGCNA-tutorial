# WGCNA-tutorial

This repository contains an R script to perform **Weighted Gene Co-expression Network Analysis (WGCNA)** using normalized RNA-seq expression data. The purpose of **WGCNA** is to identify gene modules that are co-expressed and relate them with experimental traits.  

This tutorial is based on preprocessed RNA-seq data from *Rattus norvegicus* but can be applied to other species and datasets.  

## What does it do?  

* Prepares a log-CPM expression matrix for WCGNA input  
* Identifies a soft-thresholding power for network construction
* Builds a gene co-expression network and detects gene modules
* Correlates modules with sample traits
* Visualizes module-trait relationships through heatmaps

## Requirements  

### Data 

* pre-normalized *dgeDispersion.rds* from **edgeR** (see [analysis prep pipeline](https://github.com/yasaswiveera/RNAseq-DEprep))
* *groupsVector.rds* file that assigns each sample to a group

### R packages (included in the script)  
* *WGCNA*
* *edgeR*
* *BiocManager*

## Outputs  

* WGCNAnet.rds - *WGCNA network object*
* WGCNAgenemodules.csv - *gene and module assignments*
* moduleTraitCor.csv - *module-trait correlation matrix*
* moduleTraitPval.csv - *module-trait p-value matrix*
* moduleTraitHeatmap.pdf - *heatmap of all module-trait correlations*
* moduleTraitHeatmapSignificant.pdf - *heatmap of only significant module-trait correlations (filtered by p<0.05)*
* binary heatmaps (if binary grouping is used)
