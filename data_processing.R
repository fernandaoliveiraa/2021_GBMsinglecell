
################################
# Data Processing and Analysis #
################################


# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(tidyverse)
library(readr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)


###### Load dataset ######

# Read preprocessed data from GSE131928 (https://doi.org/10.1016/j.cell.2019.06.024)
GBM_IDHwt <- read_tsv(file = "/Users/Fernanda/Desktop/R/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz")
GBM_IDHwt <- GBM_IDHwt %>% remove_rownames %>% column_to_rownames(var="GENE")

# Create the Seurat object with the non-normalized data
GBM_IDHwt <- CreateSeuratObject(counts = GBM_IDHwt, project = "GBM_IDHwt")

# Eliminate cells from pediatric patient
GBM_IDHwt <- subset(GBM_IDHwt, subset = orig.ident != "118")


###### Pre-processing and normalization ######

# MT genes were already taken off and the sequenced cells whose number of detected genes was less than half 
# or more than twice the mean number of genes detected across cells coming from the same tumor were excluded

# Apply SCtransform normalization 
GBM_IDHwt <- SCTransform(GBM_IDHwt, method = "glmGamPoi", verbose = FALSE)


###### Cell-cycle scoring and regression ######

# Segregate list of cell cycle markers (from Tirosh et al, 2015) into G2/M and S markers
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell-cycle scores
GBM_IDHwt <- CellCycleScoring(GBM_IDHwt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Regress out the difference between G2M and S phase scores and scale data 
GBM_IDHwt$CC.Difference <- GBM_IDHwt$S.Score - GBM_IDHwt$G2M.Score
GBM_IDHwt <- ScaleData(GBM_IDHwt, vars.to.regress = "CC.Difference", features = rownames(GBM_IDHwt))


###### Linear dimensional reduction (PCA) ######

# Perform PCA on the scaled data
GBM_IDHwt <- RunPCA(GBM_IDHwt, features = VariableFeatures(object = GBM_IDHwt))

# Determine the dimensionality of the dataset with ElbowPlot()
ElbowPlot(GBM_IDHwt, ndims = 50)
npcs <- 25 # choose number of pcs for downstream analysis


###### Cluster cells ######

# Construct a KNN graph based on the euclidean distance in PCA space 
GBM_IDHwt <- FindNeighbors(GBM_IDHwt, dims = 1:npcs)

# Find Clusters
resolution = 1 # set a resolution parameter
GBM_IDHwt <- FindClusters(GBM_IDHwt, resolution = resolution)

# Run UMAP
GBM_IDHwt <- RunUMAP(GBM_IDHwt, dims = 1:npcs)


###### Find cluster biomarkers ###### 

# Find markers for every cluster compared to all remaining cells, report only the positive ones
GBM_IDHwt.markers <- FindAllMarkers(GBM_IDHwt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GBM_IDHwt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


###### Assign cell type identity to clusters ###### 

# Save old identity classes
GBM_IDHwt[["old.ident"]] <- Idents(object = GBM_IDHwt)

# Rename classes
GBM_IDHwt <- RenameIdents(object = GBM_IDHwt, `0` = 'Tumor Cells', `1` = 'Tumor Cells', `2` = 'Macrophages', `3` = 'Tumor Cells', `4` = 'Macrophages',
                          `5` = 'Tumor Cells', `6` = 'Tumor Cells', `7` = 'Tumor Cells', `8` = 'Macrophages', `9` = 'Tumor Cells', `10` = 'Tumor Cells',
                          `11` = 'Macrophages', `12` = 'Macrophages', `13` = 'Tumor Cells', `14` = 'Macrophages', `15` = 'Tumor Cells', `16` = 'Macrophages',
                          `17` = 'Macrophages', `18` = 'Macrophages', `19` = 'Tumor Cells', `20` = 'Oligodendrocytes',  `21` = 'T Cells',
                          `22` = 'Tumor Cells', `23` = 'Oligodendrocytes', `24` = 'Tumor Cells', `25` = 'Tumor Cells', `26` = 'Tumor Cells', `27` = 'Tumor Cells',
                          `28` = 'Macrophages', `29` = 'Tumor Cells', `30` = 'Tumor Cells', `31` = 'Tumor Cells')




