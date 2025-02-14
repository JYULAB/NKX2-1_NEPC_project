library(copykat)
library(Seurat)

# d2
seu_d2_raw <- Read10X("../../../../../multiome/arc_ranger/d2/filtered_feature_bc_matrix", strip.suffix = TRUE)
d2 <- CreateSeuratObject(counts = seu_raw$`Gene Expression`, project = "D2", min.cells = 0, min.feature = 0)

raw_d2 <- as.matrix(d2@assays$RNA@counts)

copykat_res_d2 <- copykat(raw_d2, sam.name = "D2", cell.line="yes")
saveRDS(copykat_res_d2, "d2_copykat.RDS")


# d14
seu_d14_raw <- Read10X("../../../../../multiome/arc_ranger/d14/filtered_feature_bc_matrix", strip.suffix = TRUE)
d14 <- CreateSeuratObject(counts = seu_raw$`Gene Expression`, project = "D14", min.cells = 0, min.feature = 0)

raw_d14 <- as.matrix(d14@assays$RNA@counts)

copykat_res_d14 <- copykat(raw_d14, sam.name = "14", cell.line="yes")
saveRDS(copykat_res_d14, "d2_copykat.RDS")

# d21
seu_d21_raw <- Read10X("../../../../../multiome/arc_ranger/d21_add/filtered_feature_bc_matrix", strip.suffix = TRUE)
d21 <- CreateSeuratObject(counts = seu_raw$`Gene Expression`, project = "D2", min.cells = 0, min.feature = 0)

raw_d21 <- as.matrix(d2@assays$RNA@counts)

copykat_res_d21 <- copykat(raw_d21, sam.name = "D21", cell.line="yes")
saveRDS(copykat_res_d21, "d21_copykat.RDS")
