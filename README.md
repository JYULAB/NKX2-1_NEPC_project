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
- [**04_GEX_integration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/04_GEX_integration.R): Integration and analysis of scRNA-Seq, Fig 3C, D, E, F
- [**05_ATAC_integration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/05_ATAC_integration.R): Integration and analysis of scATAC-Seq, Fig 3H, I, J
- [**clinical_patients_intergration.R**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Fig 3D, E
- [**RNA_velocity.py**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/RNA_velocity.py), [**RNA_velocity_metadata.py**](https://github.com/JYULAB/NKX2-1_NEPC_project/blob/main/Multiome/RNA_velocity_metadata.py), Fig 3G


 
