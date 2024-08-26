# Code for NKX2-1 Paper
Citation:

This repository includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the paper.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Files needed to start analysis such as .mcool files and .hd5 files which have a large size were uploaded to GEO. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis. Expected output of each script are the indicated figures from the paper.
 
### RNA-Seq

- [**timecourse_deseq2.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/RNA-seq/timecourse_deseq2.R): FOXA2 overexpression RNA-seq Analysis, Fig 2A
- [**AR_NE_RNAseq_GSVA.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/RNA-seq/AR_NE_RNAseq_GSVA.R): AR and NR signature GSVA analysis of clinical samples, Fig 2B
- [**LuNE_DEGs.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/RNA-seq/LuNE_DEGs.R): Differential gene expression in LuNE cells, Fig 7E, 7F
- [**fpkm_pval.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/RNA-seq/fpkm_pval.R): FPKM plot of FOXA2 and NKX2-1 in human samples, Fig 5B

### Multiome

- [**01_d0_single.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/01_d0_single.R),  [**02_d14_single.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/02_d14_single.R),  [**03_d21_single.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/03_d21_single.R): Individual timepoint processing and filtering
- [**04_GEX_integration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/04_GEX_integration.R): Integration and analysis of scRNA-Seq, Fig 3C, D, E, F, Fig S3F 
- [**05_ATAC_integration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/05_ATAC_integration.R): Integration and analysis of scATAC-Seq, Fig 3H, I, J
- [**clinical_patients_intergration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Fig 3D, E
- [**RNA_velocity.py**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/RNA_velocity.py), [**RNA_velocity_metadata.py**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/RNA_velocity_metadata.py), Fig 3G
- [**Multiome_cells.csv**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/Multiome_cells.csv): cell barcodes that passed filtering

### Multiome Additional Sequencing

- [**d21_add.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/): D21 additional sequencing processing
- [**d21_add_ectopic.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/): D21 additional sequencing with inclusion of ectopic sequence processing and imputation
- [**endogenous_integration_rna.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/endogenous_integration_rna.R): Integation of RNA of original D0 and D14 with D21 additional sequencing, Fig S4B, C
- [**endogenous_integration_atac.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/endogenous_integration_atac.R): Integration of original D0 and D14 with D21 additional sequencing, Fig S4E, F
- [**meta_highlight_manual.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/meta_highlight_manual.R), [**Meta_Highlight_Plot_manual.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/Meta_Highlight_Plot_manual.R): manually tweaking plotting functions to ensure rare clones are shown on top in the UMAP, Fig S4B, E

#### copyKat
- [**01_copykat_run.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/01_copykat_run.R): Running copyKat on individual samples
- [**02_copykat_analysis.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/02_copykat_analysis.R): Filtering cells and combining individual CNA 
- [**03_draw_heatmap.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/03_draw_heatmap.R): Draw CNA heatmap, Fig S4A
- [**Multiome_additional_cells.csv**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/Multiome_additional_cells.csv): cell barcodes that passed filtering

#### epiAneufinder
- [**01_epiAneufinder_fragment_run.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/01_epiAneufinder_fragment_run.R): Running epiAneufinder on individual samples
- [**02_endogenous_epi_bins_loose.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/02_endogenous_epi_bins_loose.R): Filtering cells and combining individual CNA 
- [**03_split_subclones_color.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/03_split_subclones_color.R): Finding clones in the combined matrix and plotting heatmap, Fig S4D
- [**04_epi_loose_atac.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/04_epi_loose_atac.R): Plotting clones on UMAP and violin plots, Fig SE, F

### Hi-C

- [**Expected values.ipynb**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/Expected values.ipynb):Calculating eigenvector values by 100kb bins 
- [**AB Compartments.ipynb**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/AB Compartments.ipynb): Calculating eigenvector values by 100kb bins
- [**AB_compartment_correlation.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/AB_compartment_correlation.R): Calculating correlations between top variable compartments, Fig S1A, 2C

#### pdx_loops
Specific loops in PDXs: Fig 1A, 1B
-[**script_mustache_93.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/script_mustache_93.sh): Example script to call pairwise diffMustache
- [**script_homer_93.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/script_homer_93.sh): Example script of finding sample specific loops
- [**permutation.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/permutation.R): Permutation test for significance of CRE overlaping with specific loops, Fig 1I,J

##### nepc/crpc
- [**NEPC_loop.bedpe**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/NEPC_loop.bedpe): Specific loops file Fig 1A, 1B
- [**specific_loop.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops and finding genes linked to specific loop anchors and performing GO, Fig 1C, D, E, F
- [**script_apa**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/): Plotting APA plots Fig 1A, 1B
- [**script_homer_nepc.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/script_homer_nepc.sh): Script for finding NEPC specific loops: loops present in 2 or more samples
- [**specific_loop_cre_overlap.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/): Finding specific loops with cre overlap, Fig 1I,J
- [**script_motif.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/pdx_loops/pdx_specific_loops/nepc/script_motif.sh): Finding motifs of CREs located at specific loop anchor. For Fig 1I, J

#### lncap_loops
- [**script_mustache.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/script_mustache.sh): Script to call pairwise diffMustache between D0 and D28, Fig 2D
- [**w0_vs_w4.diffloop1.bedpe**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/w0_vs_w4.diffloop1.bedpe): Luminal and NE loops, Fig 2D 
- [**script_apa.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/script_apa.sh): Plotting APA plots Fig 2D
- [**Luminal_NE_loops.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/Luminal_NE_loops.R): Genes linked to Luminal and NE loops, S2H

##### FOXA1-NKX2-1
- [**fx_nk_liked_loops.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/FOXA1-NKX2-1/fx_nk_liked_loops.R): Finding loops with either FOXA2 or NKX2-1 at loop anchors.
- [**compare.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/FOXA1-NKX2-1/compare.sh): selecting only loops present in d28 and not d0.
- [**w4_fx_nx_no_w0_loops.bedpe**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Hi-C/lncap_loops/FOXA1-NKX2-1/w4_fx_nx_no_w0_loops.bedpe): loops used for Fig 4F

### ATAC-Seq
- [**merge_summit_by_score.sh**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/ATAC-seq/merge_summit_by_score.sh): creating consensus peaks
- [*Differentially_accessible_regions.R*](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/ATAC-seq/Differentially accessible regions.R): Finding differentially accessible regions, Fig 3A, B
- [*ATAC-seq_PCA*](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/ATAC-seq/): Including other ATAC-seq data, Fig S3C




