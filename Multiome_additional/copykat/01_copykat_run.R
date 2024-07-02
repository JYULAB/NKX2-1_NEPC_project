library(copykat)
library(Seurat)

# d2
seu_d2_raw <- Read10X_h5("../../../multiome/arc_ranger/d2/filtered_feature_bc_matrix.h5")
raw_d2 <- as.matrix(seu_d2_raw$`Gene Expression`)

d2_copykat <- copykat(raw_d2, sam.name = "D2")
saveRDS(d2_copykat, "copykat.RDS")

# d14
seu_raw <- Read10X_h5("../../../../multiome/arc_ranger/d14/filtered_feature_bc_matrix.h5")
raw <- as.matrix(seu_raw$`Gene Expression`)

copykat_res <- copykat(raw, sam.name = "D14")
saveRDS(copykat_res, "copykat_res.RDS")

# d21
seu_d21_raw <- Read10X_h5("../../../../multiome/arc_ranger/d21_add/filtered_feature_bc_matrix.h5")
raw_d21 <- as.matrix(seu_d21_raw$`Gene Expression`)

sink("R output.txt")
copykat_res <- copykat(raw_d21, sam.name = "D21_add")
saveRDS(copykat_res, "copykat_res.RDS")
