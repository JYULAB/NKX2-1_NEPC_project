# Code for NKX2-1 Paper
Citation:

This repository includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the paper.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Files needed to start analysis such as .mcool files and .hd5 files which have a large size were uploaded to GEO. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis. Expected output of each script are the indicated figures from the paper.
 
### RNA-Seq

- [**timecourse_deseq2.R**](RNA-seq/timecourse_deseq2.R): FOXA2 overexpression RNA-seq Analysis, Fig 2A
- [**AR_NE_RNAseq_GSVA.R**](RNA-seq/AR_NE_RNAseq_GSVA.R): AR and NR signature GSVA analysis of clinical samples, Fig 2B
- [**LuNE_DEGs.R**](RNA-seq/LuNE_DEGs.R): Differential gene expression in LuNE cells, Fig 7E, 7F
- [**fpkm_pval.R**](RNA-seq/fpkm_pval.R): FPKM plot of FOXA2 and NKX2-1 in human samples, Fig 5B

### Multiome

- [**01_d0_single.R**](Multiome/01_d0_single.R),  [**02_d14_single.R**](Multiome/02_d14_single.R),  [**03_d21_single.R**](Multiome/03_d21_single.R): Individual timepoint processing and filtering
- [**04_GEX_integration.R**](Multiome/04_GEX_integration.R): Integration and analysis of scRNA-Seq, Fig 3C, D, E, F, Fig S3F 
- [**05_ATAC_integration.R**](Multiome/05_ATAC_integration.R): Integration and analysis of scATAC-Seq, Fig 3H, I, J
- [**clinical_patients_intergration.R**](Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Fig 3D, E
- [**RNA_velocity.py**](Multiome/RNA_velocity.py), [**RNA_velocity_metadata.py**](Multiome/RNA_velocity_metadata.py), Fig 3G
- [**Multiome_cells.csv**](Multiome/Multiome_cells.csv): cell barcodes that passed filtering

### Multiome Additional Sequencing

- [**d21_add.R**](Multiome_additional/): D21 additional sequencing processing
- [**d21_add_ectopic.R**](Multiome_additional/): D21 additional sequencing with inclusion of ectopic sequence processing and imputation
- [**endogenous_integration_rna.R**](Multiome_additional/endogenous_integration_rna.R): Integation of RNA of original D0 and D14 with D21 additional sequencing, Fig S4B, C
- [**endogenous_integration_atac.R**](Multiome_additional/endogenous_integration_atac.R): Integration of original D0 and D14 with D21 additional sequencing, Fig S4E, F
- [**meta_highlight_manual.R**](Multiome_additional/meta_highlight_manual.R), [**Meta_Highlight_Plot_manual.R**](Multiome_additional/Meta_Highlight_Plot_manual.R): manually tweaking plotting functions to ensure rare clones are shown on top in the UMAP, Fig S4B, E

#### copyKat
- [**01_copykat_run.R**](Multiome_additional/01_copykat_run.R): Running copyKat on individual samples
- [**02_copykat_analysis.R**](Multiome_additional/02_copykat_analysis.R): Filtering cells and combining individual CNA 
- [**03_draw_heatmap.R**](Multiome_additional/03_draw_heatmap.R): Draw CNA heatmap, Fig S4A
- [**Multiome_additional_cells.csv**](Multiome_additional/Multiome_additional_cells.csv): cell barcodes that passed filtering

#### epiAneufinder
- [**01_epiAneufinder_fragment_run.R**](Multiome_additional/01_epiAneufinder_fragment_run.R): Running epiAneufinder on individual samples
- [**02_endogenous_epi_bins_loose.R**](Multiome_additional/02_endogenous_epi_bins_loose.R): Filtering cells and combining individual CNA 
- [**03_split_subclones_color.R**](Multiome_additional/03_split_subclones_color.R): Finding clones in the combined matrix and plotting heatmap, Fig S4D
- [**04_epi_loose_atac.R**](Multiome_additional/04_epi_loose_atac.R): Plotting clones on UMAP and violin plots, Fig SE, F

### Hi-C  
- [**Expected_values.ipynb**](Hi-C/Expected_values.ipynb):Calculating eigenvector values by 100kb bins 
- [**AB_Compartments.ipynb**](Hi-C/AB_Compartments.ipynb): Calculating eigenvector values by 100kb bins
- [**AB_compartment_correlation.R**](Hi-C/AB_compartment_correlation.R): Calculating correlations between top variable compartments, Fig S1A, 2C

#### pdx_loops
Specific loops in PDXs: Fig 1A, 1B  
- [**script_mustache_93.sh**](Hi-C/pdx_loops/script_mustache_93.sh): Example script to call pairwise diffMustache
- [**script_homer_93.sh**](Hi-C/pdx_loops/pdx_specific_loops/script_homer_93.sh): Example script of finding sample specific loops
- [**permutation.R**](Hi-C/pdx_loops/pdx_specific_loops/permutation.R): Permutation test for significance of CRE overlaping with specific loops, Fig 1I,J

##### nepc/crpc
- [**NEPC_loop.bedpe**](Hi-C/pdx_loops/pdx_specific_loops/nepc/NEPC_loop.bedpe): Specific loops file Fig 1A, 1B
- [**specific_loop.R**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops and finding genes linked to specific loop anchors and performing GO, Fig 1C, D, E, F
- [**script_apa**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Plotting APA plots Fig 1A, 1B
- [**script_homer_nepc.sh**](Hi-C/pdx_loops/pdx_specific_loops/nepc/script_homer_nepc.sh): Script for finding NEPC specific loops: loops present in 2 or more samples
- [**specific_loop_cre_overlap.R**](Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops with cre overlap, Fig 1I,J
- [**script_motif.sh**](Hi-C/pdx_loops/pdx_specific_loops/nepc/script_motif.sh): Finding motifs of CREs located at specific loop anchor. For Fig 1I, J

#### lncap_loops
- [**script_mustache.sh**](Hi-C/lncap_loops/script_mustache.sh): Script to call pairwise diffMustache between D0 and D28, Fig 2D
- [**w0_vs_w4.diffloop1.bedpe**](Hi-C/lncap_loops/w0_vs_w4.diffloop1.bedpe): Luminal and NE loops, Fig 2D 
- [**script_apa.sh**](Hi-C/lncap_loops/script_apa.sh): Plotting APA plots Fig 2D
- [**Luminal_NE_loops.R**](Hi-C/lncap_loops/Luminal_NE_loops.R): Genes linked to Luminal and NE loops, S2H

##### FOXA1-NKX2-1
- [**fx_nk_liked_loops.R**](Hi-C/lncap_loops/FOXA1-NKX2-1/fx_nk_liked_loops.R): Finding loops with either FOXA2 or NKX2-1 at loop anchors
- [**compare.sh**](Hi-C/lncap_loops/FOXA1-NKX2-1/compare.sh): selecting only loops present in d28 and not d0
- [**w4_fx_nx_no_w0_loops.bedpe**](Hi-C/lncap_loops/FOXA1-NKX2-1/w4_fx_nx_no_w0_loops.bedpe): loops used for Fig 4F

### ATAC-seq
- [**merge_summit_by_score.sh**](ATAC-seq/merge_summit_by_score.sh): creating consensus peaks  
- [**Differentially_accessible_regions.R**](ATAC-seq/Differentially_accessible_regions.R): Finding differentially accessible regions, Fig 3A, B
- [**ATAC-seq_PCA.R**](ATAC-seq/ATAC-seq_PCA.R): Including other ATAC-seq data, Fig S3C

### ChIP-seq
- [**command.sh**](ChIP-seq/command.sh): creating heatmaps, Fig 4D, 4E, 4H, 5E, 7D, S5E, S7E, S9A, S9B (same tool is used for ATAC-seq heatmaps).

### miscellaneous_code
- [**gene_correlation_for_all_tables_4ppt_ggplot.R**](miscellaneous_code/gene_correlation_for_all_tables_4ppt_ggplot.R): R code to draw scatter plots with prostate cancer patient data, Fig S6A, S7B  
- [**command.sh**](miscellaneous_code/command.sh): generating super enhancer plots, Fig S9C, S9E
- [**TF_volcano.R**](miscellaneous_code/TF_volcano.R): creating TF volcano plots, Fig 5A

### Nanopore_long_read
- [**DNA_methylation_boxplot.R**](Nanopore_long_read/DNA_methylation_boxplot.R): R code to draw patient DNA methylation level box plots  Fig S6C
- [**IGV_track_command.sh**](Nanopore_long_read/IGV_track_command.sh): generating IGV bigwig track files, Fig 6A, 6E, S6B
- [**command.sh**](Nanopore_long_read/command.sh): generating mA and mCpG heatmaps and profiles, Fig 6B, S8A
- [**mA_mCpG_corr_plot.R**](Nanopore_long_read/mA_mCpG_corr_plot.R): creating mA and mCpG correlation plots, Fig 6C, 6G
- [**methylation_PCA.R**](Nanopore_long_read/methylation_PCA.R): creating methylation profile PCA plots, Fig 6D, S8E
- [**PDX_methylation_command.sh**](Nanopore_long_read/PDX_methylation_command.sh): generating mA and mCpG of PDX profiles, Fig 6F, S8B, S8C
