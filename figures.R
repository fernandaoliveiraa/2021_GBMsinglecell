
################################
#            Figures           #
################################


# Load packages
library(RColorBrewer)
library(ggpubr)
library(scales)
library(ggridges)


###### Figure 1 ######


# Figure 1A

# DimPlot
DimPlot(GBM_IDHwt, label = FALSE, pt.size = 0.85) + border("black", size = 1) + 
  theme(legend.position = c(0.7, 0.9),
        axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(colour="black", size=15, margin = margin(t = 5)))

# List of marker genes
marker_genes <- c('MBP', 'TF', 'PLP1', 'MAG', 'MOG', 'CLDN11','CD14', 'AIF1', 'FCER1G', 'FCGR3A', 'TYROBP', 
                  'CSF1R','CD2', 'CD3D', 'CD3E', 'CD3G','SOX2', 'PARP1')

# DotPlot
DotPlot(GBM_IDHwt, features = marker_genes, cols = rev(brewer.pal(n = 11, name = "Spectral")), dot.scale = 20) + border("black", size = 0.8) + 
  theme(axis.line = element_line(size = 0.8), axis.title = element_blank(), axis.text.y = element_text(size = 25), 
        axis.text.x = element_text(size = 20), axis.ticks = element_blank()) + RotatedAxis() + coord_flip() + NoLegend()


# Figure 1B

DimPlot(GBM_IDHwt, label = FALSE, pt.size = 0.85, group.by = 'orig.ident') + border("black", size = 1) +
  theme(legend.position = c(0.88, 0.8),
        axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(colour="black", size=15, margin = margin(t = 5)),
        plot.title = element_blank())


# Figure 1C

# Er stress signature gene list
erstress_gene_signature <- c('DDIT3','CHAC1','HERPUD1','DNAJB9','CTH','SESN2','XBP1','EIF2AK3','ATF3','BHLHE40', 'NFIL3','PIGA','TSPYL2', 'MXD1','MIS12','RAB39B','NR4A1','PELI1','PTGS2','NEDD9','EGR3','RCAN1','PDE4D','LYSMD3','EGR2','BTG2','TRIB3','SLC10A7','WIPI1','TUBE1','SERP1','PGM3','CCDC174','NRBF2','ATP6V1C2','SNX22','CDK2AP2')

# Create a metafeature with the average expression of the gene list
GBM_IDHwt <- MetaFeature(object = GBM_IDHwt, features = erstress_gene_signature, meta.name = "erstress_gene_signature")

# FeaturePlot
FeaturePlot(object = GBM_IDHwt,  pt.size = 0.85, features = "erstress_gene_signature", label = FALSE, min.cutoff = "q10", max.cutoff = "q90") +
  border("black", size = 1) + 
  theme(axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        plot.title = element_blank())


# Figure 1D

RidgePlot(GBM_IDHwt, features = "erstress_gene_signature") + 
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank(),
        legend.position = 'bottom', legend.text = element_text(colour="black", size=15), legend.spacing.x = unit(10, 'pt'))


# Figure 1E

FeaturePlot(object = GBM_IDHwt,  pt.size = 0.85, features = c('XBP1', 'HERPUD1', 'P4HB', 'HSPA5'), 
            label = FALSE, min.cutoff = "q10", max.cutoff = "q90") & 
  border("black", size = 0.8) & 
  theme(axis.line = element_line(size = 0.8),
        axis.title = element_text(size = 15)) 


# Figure 1F

RidgePlot(GBM_IDHwt, features = c('XBP1', 'HERPUD1', 'P4HB', 'HSPA5'), ncol = 2) &
  theme_ridges(grid = TRUE, line_size = 0.7) & border("black", size = 0.8) & 
  theme(axis.text.y = element_blank(), axis.title = element_blank()) &
  NoLegend()


###### Figure S1 ######


# Figure S1A

# Tumor Cells
RidgePlot(GBM_IDHwt, features = "erstress_gene_signature", group.by = 'orig.ident', idents = 'Tumor Cells') + 
  xlim(-0.00005, 0.0004) +
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank(),
        legend.position = 'right', legend.text = element_text(colour="black", size=15), legend.spacing.x = unit(10, 'pt'))

# Macrophages
RidgePlot(GBM_IDHwt, features = "erstress_gene_signature", group.by = 'orig.ident', idents = 'Macrophages') +
  xlim(-0.00005, 0.0004) +
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank(),
        legend.position = 'right', legend.text = element_text(colour="black", size=15), legend.spacing.x = unit(10, 'pt')) 


# Figure S1B

# Create a Seurat object with only cells from one patient
p105 <- subset(GBM_IDHwt, subset = orig.ident == '105')

# FeaturePlot
FeaturePlot(object = p105,  pt.size = 0.85, features = "erstress_gene_signature", label = FALSE, 
            min.cutoff = 0.00001, max.cutoff = 0.00016) +
  border("black", size = 1) + 
  theme(axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        plot.title = element_blank())

# RidgePlot
RidgePlot(p105, features = "erstress_gene_signature", idents = c('Tumor Cells', 'Macrophages')) +  xlim(-0.00005, 0.0005) +
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank(),
        legend.position = 'top', legend.text = element_text(colour="black", size=15), legend.spacing.x = unit(10, 'pt'))


# Figure S1C

# Create a Seurat object with only cells from one patient
p115 <- subset(GBM_IDHwt, subset = orig.ident == '115')

# FeaturePlot
FeaturePlot(object = p115,  pt.size = 0.85, features = "erstress_gene_signature", label = FALSE, 
            min.cutoff = 0.00001, max.cutoff = 0.00016) +
  border("black", size = 1) + 
  theme(axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        plot.title = element_blank())

# RidgePlot
RidgePlot(p115, features = "erstress_gene_signature", idents = c('Tumor Cells', 'Macrophages')) +  xlim(-0.00005, 0.0005) +
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank()) + NoLegend()


# Figure S1D

# Create a Seurat object with only cells from one patient
p124 <- subset(GBM_IDHwt, subset = orig.ident == '115')

# FeaturePlot
FeaturePlot(object = p124,  pt.size = 0.85, features = "erstress_gene_signature", label = FALSE, 
            min.cutoff = 0.00001, max.cutoff = 0.00016) +
  border("black", size = 1) + 
  theme(axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        plot.title = element_blank())

# RidgePlot
RidgePlot(p124, features = "erstress_gene_signature", idents = c('Tumor Cells', 'Macrophages')) +  xlim(-0.00005, 0.0005) +
  theme_ridges(grid = TRUE, line_size = 0.7) + 
  theme(axis.title = element_blank(), plot.title = element_blank(), axis.text.y = element_blank()) + NoLegend()


###### Figure S3 ######


# Figure S3C

# Define a pseudocount to add to the data in order to avoid Inf and -Inf values
pseudocount <- 0.1

# Compute ratio of expression values
gene.exp.ratio <- (GBM_IDHwt@assays[["SCT"]]@data["XBP1",] + pseudocount) / (GBM_IDHwt@assays[["SCT"]]@data["MAPK8",] + pseudocount)

# Add ratio of expression values in meta.data
if (all(names(x = gene.exp.ratio) == rownames(x = GBM_IDHwt@meta.data))) {
  cat("Cell names order match in 'gene.exp.ratio' and 'object@meta.data':\n", 
      "adding ratio of expression values in 'object@meta.data$gene.exp.ratio'")
  GBM_IDHwt@meta.data$gene.exp.ratio <- gene.exp.ratio
}

# FeaturePlot
FeaturePlot(object = GBM_IDHwt,  pt.size = 0.85, features = "gene.exp.ratio", label = FALSE) +
  border("black", size = 1) + 
  theme(axis.line = element_line(size = 1),
        axis.title = element_text(size = 20),
        plot.title = element_blank())

