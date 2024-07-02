# Code for NKX2-1 Paper
Citation:

This repositopry includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the manuscript.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Files needed to start analysis such as .mcool files and .hd5 files which have a large size were uploaded to GEO. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis.
Expected output of each script are the indicated figures from the paper.
 
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
- [**endogenous_integration_rna.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/_additional/endogenous_integration_rna.R): Integation of RNA of original D0 and D14 with D21 additional sequencing, Fig S4B, C
- [**endogenous_integration_atac.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/_additional/endogenous_integration_atac.R): Integration of original D0 and D14 with D21 additional sequencing, Fig S4E, F
- [**meta_highlight_manual.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/meta_highlight_manual.R), [**Meta_Highlight_Plot_manual.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/Meta_Highlight_Plot_manual.R): manually tweaking plotting functions to ensure rare clones are shown on top in the UMAP, Fig S4B, E

#### copyKat
- [**01_copykat_run.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/01_copykat_run.R): Running copyKat on individual samples
- [**02_copykat_analysis.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/02_copykat_analysis.R): Filtering cells and combining individual CNA 
- [**03_draw_heatmap.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/03_draw_heatmap.R): Draw CNA heatmap, Fig S4A

#### epiAneufinder
- [**01_epiAneufinder_fragment_run.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/01_epiAneufinder_fragment_run.R): Running epiAneufinder on individual samples
- [**02_endogenous_epi_bins_loose.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/02_endogenous_epi_bins_loose.R): Filtering cells and combining individual CNA 
- [**03_split_subclones_color.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/03_split_subclones_color.R): Finding clones in the combined matrix and plotting heatmap, Fig S4D
- [**04_epi_loose_atac.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome_additional/04_epi_loose_atac.R): Plotting clones on UMAP and violin plots, Fig SE, F

 
