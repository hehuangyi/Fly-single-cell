This repository contains the workflow and codes of the major analyse in the manuscript, titled "New genes play a prominent role in evolution of new sperm classes of Drosophila". 

# Figure plot

This folder contains scripts used to generate the main figures presented in the manuscript.

# New gene process

This folder contains the pipeline for identifying different types of new genes and input files. 

- `01.duplication_translocation_detection.R`    Define Dpse duplication and translocation events
    
- `02.singlecopy_detection.R`    Identify and extract Dpse single-copy orthologs

- `03.denovo_gene_blast.sh`    Runs BLASTP and TBLASTN of Dpse proteins across species and extracts best hits

-  `04.denovo_gene_genewise_spaln.sh`    Performs genome alignment, preprocessing, and cross-species homology annotation (Genewise/Spaln) to identify and validate de novo genes

-  `05.denovo_filter_pipeline.py`    Identify high-confidence de novo genes by intersecting BLAST- and Spaln-based filters

# RNA seq process
This folder contains scripts for RNA-seq data processing and expression analysis
- `RNAseq_process.sh`    Pipeline for processing RNA-seq data
- `RPKM_and_testis_specificity_calculate.R` Calculatie gene expression levels (RPKM) and assessing testis specificity of genes
    
# scRNA-seq process
-  `run_scRNA-seq.sh`    Run CellRanger to get count matrix
- `01.preprocess.R`    QC filtering, normalization, and clustering for single-cell RNA-seq data
- `02.doubletfinder.R`    Define doublet using DoubletFinder
- `03.cross_species_integration.R`    Integrate cells from different stages and species
- `04.Dpse_cell_integration.R`    Integrate Dpse cells from different stages

# ChIP-seq process
- `run_chipseq.sh`    Processes alignment, normalization of ChIP-seq data and generate gene-level IP/Input files.

# data
This folder contains input files used by the analysis scripts.

If any questions, please contract Huangyi He: hehuangyi@zju.edu.cn or Prof. Qi Zhou: zhouqi1982@zju.edu.cn.
