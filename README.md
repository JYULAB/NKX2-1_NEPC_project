# Code for NKX2-1 Paper
Citation:

This repository includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the paper.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Files needed to start analysis such as .mcool files and .hd5 files which have a large size were uploaded to GEO. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis. Expected output of each script are the indicated figures from the paper.
 
### RNA-Seq

- [**timecourse_deseq2.R**](RNA-seq/timecourse_deseq2.R): FOXA2 overexpression RNA-seq Analysis, Fig. 2a
- [**AR_NE_RNAseq_GSVA.R**](RNA-seq/AR_NE_RNAseq_GSVA.R): AR and NR signature GSVA analysis of clinical samples, Fig. 2c
- [**LuNE_DEGs.R**](RNA-seq/LuNE_DEGs.R): Differential gene expression in LuNE cells, Fig. 7e, 7f
- [**fpkm_pval.R**](RNA-seq/fpkm_pval.R): FPKM plot of FOXA2 and NKX2-1 in human samples, Fig. 5b

### Multiome

- [**01_d0_single.R**](Multiome/01_d0_single.R),  [**02_d14_single.R**](Multiome/02_d14_single.R),  [**03_d21_single.R**](Multiome/03_d21_single.R): Individual timepoint processing and filtering
- [**04_GEX_integration.R**](Multiome/04_GEX_integration.R): Integration and analysis of scRNA-Seq, Fig. 3c, 3d, 3e, 3f, Extended Data Fig. 3f 
- [**05_ATAC_integration.R**](Multiome/05_ATAC_integration.R): Integration and analysis of scATAC-Seq, Fig. 3h, 3i, 3j
- [**clinical_patients_intergration.R**](Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Extended Data Fig. 3d, 3e
- [**RNA_velocity.py**](Multiome/RNA_velocity.py), [**RNA_velocity_metadata.py**](Multiome/RNA_velocity_metadata.py), Fig. 3g
- [**Multiome_cells.csv**](Multiome/Multiome_cells.csv): cell barcodes that passed filtering

### Multiome Additional Sequencing

- [**d21_add.R**](Multiome_additional/d21_add.R): D21 additional sequencing processing
- [**d21_add_ectopic.R**](Multiome_additional/d21_add_ectopic.R): D21 additional sequencing with inclusion of ectopic sequence processing and imputation
- [**meta_highlight_manual.R**](Multiome_additional/meta_highlight_manual.R), [**Meta_Highlight_Plot_manual.R**](Multiome_additional/Meta_Highlight_Plot_manual.R): manually tweaking plotting functions to ensure rare clones are shown on top in the UMAP, Extended Data Fig. 4b, 4e

#### copyKat
- [**01_copykat_run.R**](Multiome_additional/copykat/01_copykat_run.R): Running copyKat on individual samples
- [**02_copykat_analysis.R**](Multiome_additional/copykat/02_copykat_analysis.R): Filtering cells and combining individual CNA 
- [**03_draw_heatmap.R**](Multiome_additional/copykat/03_draw_heatmap.R): Draw CNA heatmap, Extended Data Fig. 4a
- [**04_rna_integration_cna_filtered.R**](Multiome_additional/copykat/04_rna_integration_cna_filtered.R): Integration of individual timepoints and plotting, Extended Data Fig. 4b, 4c
- [**Multiome_additional_cells.csv**](Multiome_additional/copykat/copy_cells.csv): cell barcodes that passed filtering
  
#### scNanoRNA-seq
- [**d0-d28_scnanornaseq.R**](scNanoRNA-seq/d0-d28_scnanornaseq.R): Script for scNanoRNA-seq analysis, Extended Data Fig. 4d-h

### Hi-C  
- [**Expected_values.ipynb**](Hi-C/Expected_values.ipynb):Calculating eigenvector values by 100kb bins 
- [**AB_Compartments.ipynb**](Hi-C/AB_Compartments.ipynb): Calculating eigenvector values by 100kb bins
- [**AB_compartment_correlation.R**](Hi-C/AB_compartment_correlation.R): Calculating correlations between top variable compartments, Extended Data Fig. 1a, 1g

#### pdx_loops
Specific loops in PDXs: Fig. 1a, 1b  
- [**script_mustache_93.sh**](Hi-C/pdx_loops/script_mustache_93.sh): Example script to call pairwise diffMustache
- [**script_homer_93.sh**](Hi-C/pdx_loops/pdx_specific_loops/script_homer_93.sh): Example script of finding sample specific loops
- [**permutation.R**](Hi-C/pdx_loops/pdx_specific_loops/permutation.R): Permutation test for significance of CRE overlaping with specific loops, Fig. 1i,j

##### nepc/crpc
- [**NEPC_loop.bedpe**](Hi-C/pdx_loops/pdx_specific_loops/nepc/NEPC_loop.bedpe): Specific loops file Fig. 1a, 1b
- [**specific_loop.R**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops and finding genes linked to specific loop anchors and performing GO, Fig. 1c, 1d, 1e, 1f
- [**script_apa**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Plotting APA plots Fig. 1a, 1b
- [**script_homer_nepc.sh**](Hi-C/pdx_loops/pdx_specific_loops/nepc/script_homer_nepc.sh): Script for finding NEPC specific loops: loops present in 2 or more samples
- [**specific_loop_cre_overlap.R**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops with cre overlap, Fig. 1i,1j
- [**script_motif.sh**](Hi-C/pdx_loops/pdx_specific_loops/nepc/script_motif.sh): Finding motifs of CREs located at specific loop anchor. For Fig. 1i,1j

#### lncap_loops
- [**script_mustache.sh**](Hi-C/lncap_loops/script_mustache.sh): Script to call pairwise diffMustache between D0 and D28, Fig. 2d
- [**w0_vs_w4.diffloop1.bedpe**](Hi-C/lncap_loops/w0_vs_w4.diffloop1.bedpe): Luminal and NE loops, Fig. 2d 
- [**script_apa.sh**](Hi-C/lncap_loops/script_apa.sh): Plotting APA plots Fig. 2d
- [**Luminal_NE_loops.R**](Hi-C/lncap_loops/Luminal_NE_loops.R): Genes linked to Luminal and NE loops, Supplementary Data Fig. 1a

##### FOXA1-NKX2-1
- [**fx_nk_liked_loops.R**](Hi-C/lncap_loops/FOXA1-NKX2-1/fx_nk_liked_loops.R): Finding loops with either FOXA2 or NKX2-1 at loop anchors
- [**compare.sh**](Hi-C/lncap_loops/FOXA1-NKX2-1/compare.sh): selecting only loops present in d28 and not d0
- [**w4_fx_nx_no_w0_loops.bedpe**](Hi-C/lncap_loops/FOXA1-NKX2-1/w4_fx_nx_no_w0_loops.bedpe): loops used for Fig. 4f

### ATAC-seq
- [**merge_summit_by_score.sh**](ATAC-seq/merge_summit_by_score.sh): creating consensus peaks  
- [**Differentially_accessible_regions.R**](ATAC-seq/Differentially_accessible_regions.R): Finding differentially accessible regions, Fig. 3a, 3b
- [**ATAC-seq_PCA.R**](ATAC-seq/ATAC-seq_PCA.R): Including LuCaP and other ATAC-seq data, Extended Data Fig. 3c

### ChIP-seq
- [**command.sh**](ChIP-seq/command.sh): creating heatmaps, Fig.4d, 4e, 4h, 5e, 7d, Extended Data Fig.5e, 7e, 9a, 9b (same tool is used for ATAC-seq heatmaps).

### miscellaneous_code
- [**gene_correlation_for_all_tables_4ppt_ggplot.R**](miscellaneous_code/gene_correlation_for_all_tables_4ppt_ggplot.R): R code to draw scatter plots with prostate cancer patient data, Extended Data Fig.6a, 7b
- [**command.sh**](miscellaneous_code/command.sh): generating super enhancer plots, Extended Data Fig.9c, 9e
- [**TF_volcano.R**](miscellaneous_code/TF_volcano.R): creating TF volcano plots, Fig.5a

### Nanopore_long_read
- [**DNA_methylation_boxplot.R**](Nanopore_long_read/DNA_methylation_boxplot.R): R code to draw patient DNA methylation level box plots, Extended Data Fig.6c
- [**IGV_track_command.sh**](Nanopore_long_read/IGV_track_command.sh): generating IGV bigwig track files, Fig.6a, 6e, Extended Data Fig.6b
- [**command.sh**](Nanopore_long_read/command.sh): generating mA and mCpG heatmaps and profiles, Fig.6b, Extended Data Fig.8a
- [**mA_mCpG_corr_plot.R**](Nanopore_long_read/mA_mCpG_corr_plot.R): creating mA and mCpG correlation plots, Fig.6c, 6g
- [**methylation_PCA.R**](Nanopore_long_read/methylation_PCA.R): creating methylation profile PCA plots, Fig.6d, Extended Data Fig.8c
- [**PDX_methylation_command.sh**](Nanopore_long_read/PDX_methylation_command.sh): generating mA and mCpG of PDX profiles, Fig.6f, Extended Data Fig.8d, 8e
