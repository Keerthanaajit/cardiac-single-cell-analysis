# PACKAGES
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(scuttle)
library(DropletUtils)
library(scDblFinder)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(tibble)
library(ggrepel)

# DATA INPUT
# The sequencing data used in this study are not publicly
# available. Users wishing to reproduce the analysis should
# provide their own Parse Biosciences gene-expression matrices
# and specify the appropriate paths below.

data_dir <- "path/to/data"
# Read samples and create Seurat objects
# Sample identifiers have been replaced with generic labels.

sample_01 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_01", "DGE_unfiltered")
)
sample_01 <- CreateSeuratObject(counts = sample_01, project = "sample_01")
sample_01$sex <- "Female"
sample_01$treatment <- "Ang"
sample_01$time_point <- "10d"


sample_02 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_02", "DGE_unfiltered")
)
sample_02 <- CreateSeuratObject(counts = sample_02, project = "sample_02")
sample_02$sex <- "Male"
sample_02$treatment <- "Ang"
sample_02$time_point <- "14d"


sample_03 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_03", "DGE_unfiltered")
)
sample_03 <- CreateSeuratObject(counts = sample_03, project = "sample_03")
sample_03$sex <- "Male"
sample_03$treatment <- "Ang"
sample_03$time_point <- "14d"


sample_04 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_04", "DGE_unfiltered")
)
sample_04 <- CreateSeuratObject(counts = sample_04, project = "sample_04")
sample_04$sex <- "Male"
sample_04$treatment <- "Sal"
sample_04$time_point <- "14d"


sample_05 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_05", "DGE_unfiltered")
)
sample_05 <- CreateSeuratObject(counts = sample_05, project = "sample_05")
sample_05$sex <- "Male"
sample_05$treatment <- "Sal"
sample_05$time_point <- "10d"


sample_06 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_06", "DGE_unfiltered")
)
sample_06 <- CreateSeuratObject(counts = sample_06, project = "sample_06")
sample_06$sex <- "Female"
sample_06$treatment <- "Ang"
sample_06$time_point <- "10d"


sample_07 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_07", "DGE_unfiltered")
)
sample_07 <- CreateSeuratObject(counts = sample_07, project = "sample_07")
sample_07$sex <- "Male"
sample_07$treatment <- "Sal"
sample_07$time_point <- "14d"


sample_08 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_08", "DGE_unfiltered")
)
sample_08 <- CreateSeuratObject(counts = sample_08, project = "sample_08")
sample_08$sex <- "Female"
sample_08$treatment <- "Ang"
sample_08$time_point <- "14d"


sample_09 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_09", "DGE_unfiltered")
)
sample_09 <- CreateSeuratObject(counts = sample_09, project = "sample_09")
sample_09$sex <- "Male"
sample_09$treatment <- "Ang"
sample_09$time_point <- "10d"


sample_10 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_10", "DGE_unfiltered")
)
sample_10 <- CreateSeuratObject(counts = sample_10, project = "sample_10")
sample_10$sex <- "Male"
sample_10$treatment <- "Ang"
sample_10$time_point <- "10d"


sample_11 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_11", "DGE_unfiltered")
)
sample_11 <- CreateSeuratObject(counts = sample_11, project = "sample_11")
sample_11$sex <- "Male"
sample_11$treatment <- "Sal"
sample_11$time_point <- "10d"


sample_12 <- ReadParseBio(
  data.dir = file.path(data_dir, "sample_12", "DGE_unfiltered")
)
sample_12 <- CreateSeuratObject(counts = sample_12, project = "sample_12")
sample_12$sex <- "Female"
sample_12$treatment <- "Ang"
sample_12$time_point <- "10d"


# MERGE SAMPLES
Ang1 <- merge(sample_01,y = list(sample_02, sample_03, sample_04, sample_05, sample_06, sample_07, sample_08, sample_09, sample_10, sample_11, sample_12), project = "Ang")

# Inspect number for genes and cells before the filtering 
cat("After merge:", "Genes =", nrow(Ang1), "Cells =", ncol(Ang1), "\n")

# Inspect number of cells in treatment groups 
table(Ang1$treatment)

# EmptyDrops filtering
# Identify cell-containing droplets from the raw count matrix and retain droplets passing the specified FDR threshold.
# Join layers before converting
Ang1[["RNA"]] <- JoinLayers(Ang1[["RNA"]])

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(Ang1)

# EmptyDrops cell calling
set.seed(123)
e.out <- emptyDrops(counts(sce),by.rank = 25000,lower = 35)

# Keep cells with FDR <= 0.001
sce <- sce[, which(e.out$FDR <= 0.001)]

# Inspect the number of genes and cells after EmptryDrops filtering 
cat("After EmptyDrops:",
    "Genes =", nrow(sce),
    "Cells =", ncol(sce), "\n")

# UMI filtering
# Calculate per-cell QC metrics and remove cells with very low total UMI counts.
qc <- perCellQCMetrics(sce)
# Remove cells with <400 UMIs
qc.lib <- qc$sum < 400
discard <- qc.lib 
cat("Low UMI cells:", sum(qc.lib), "\n")
cat("Total removed:", sum(discard), "\n")
sce <- sce[, !discard]

# Inspect genes and cells after UMI filtering 
cat("After UMI filtering:",
    "Genes =", nrow(sce),
    "Cells =", ncol(sce), "\n")

# removing genes detected in fewer than 3 cells
# Filter very sparsely detected genes to reduce noise and unnecessary features in downstream analyses.
Ang1 <- as.Seurat(sce,counts = "counts",data = NULL)
counts_matrix <- GetAssayData(Ang1,assay = "RNA",layer = "counts")
genes_keep <- Matrix::rowSums(counts_matrix > 0) >= 3
Ang1 <- subset(Ang1,features = rownames(counts_matrix)[genes_keep])

# Recalculate QC metadata after gene filtering
counts_matrix <- GetAssayData(Ang1,assay = "RNA",layer = "counts")
Ang1$nCount_RNA <- Matrix::colSums(counts_matrix)
Ang1$nFeature_RNA <- Matrix::colSums(counts_matrix > 0)

cat(
  "After converting to Seurat and filtering rare genes:",
  "Genes =", nrow(Ang1),
  "Cells =", ncol(Ang1),
  "\n"
)
# Calculate the percentage of mitochondrial transcripts per cell as an additional QC metric.
Ang1[["percent.mt"]] <- PercentageFeatureSet(Ang1,pattern = "^mt-")

# Visualize QC metrics as a violin plot 
qc_sample <- VlnPlot(Ang1,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident",ncol = 3,pt.size = 0)
print(qc_sample)

plot1 <- FeatureScatter(Ang1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Ang1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
# Retain cells with 201–4,999 detected genes and less than 5% mitochondrial RNA.
Ang1.filtered <- subset(Ang1,subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
# Inspect cells after QC filtering 
cat("After Seurat QC filtering:",
    "Genes =", nrow(Ang1.filtered),
    "Cells =", ncol(Ang1.filtered), "\n")

ncol(Ang1.filtered)
table(Ang1.filtered$treatment)
table(Ang1.filtered$orig.ident)

# Normalizing the data 
Ang1.filtered <- NormalizeData (Ang1.filtered, normalization.method = 'LogNormalize', scale.factor = 10000)

# Identifying highly variable genes (feature selection)
Ang1.filtered <- FindVariableFeatures(Ang1.filtered, selection.method ='vst', nfeatures = 2000)

# Top 10 highly variabele genes 
# Select the 2,000 genes showing the strongest standardized expression variability for dimensional reduction.
top10 <- head(VariableFeatures(Ang1.filtered),10)
# Plot variable features with and without labels 
plot1 <- VariableFeaturePlot(Ang1.filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1 + plot2)

# Scaling the data (all genes) 
all.genes <- rownames(Ang1.filtered)
Ang1.filtered <- ScaleData(Ang1.filtered, features = all.genes)

# linear diemntional reduction 
Ang1.filtered <- RunPCA(Ang1.filtered, features = VariableFeatures(object = Ang1.filtered))

# Examine and visualize PCA results a few different ways 
print(Ang1.filtered[['pca']], dims = 1:5, nfeatures = 5)

VizDimLoadings(Ang1.filtered, dims = 1:2, reduction ='pca')

DimPlot(Ang1.filtered, reduction = 'pca') + NoLegend()

# Display the PCA heatmaps                            
DimHeatmap(Ang1.filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Ang1.filtered, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot 
ElbowPlot(Ang1.filtered, ndims = 30)

# Cluster the cells 
# Build nearest-neighbor graphs and compare clustering using the first 10 versus first 15 principal components. 
# PC1-10 clustering
Ang1.PC10 <- FindNeighbors(Ang1.filtered, dims = 1:10)
Ang1.PC10 <- FindClusters(Ang1.PC10, resolution = 0.3)

# PC1-15 clustering
Ang1.PC15 <- FindNeighbors(Ang1.filtered, dims = 1:15)
Ang1.PC15 <- FindClusters(Ang1.PC15, resolution = 0.3)

# Look at cluster IDs of the first 5 cells
head(Idents(Ang1.PC10), 5)
head(Idents(Ang1.PC15), 5)

# UMAP Clustering 
set.seed(123)
Ang1.PC10 <- RunUMAP(Ang1.PC10, dims = 1:10)
set.seed(123)
Ang1.PC15 <- RunUMAP(Ang1.PC15, dims = 1:15)

DimPlot(Ang1.PC10, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(Ang1.PC15, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

length(levels(Idents(Ang1.PC10)))
length(levels(Idents(Ang1.PC15)))

table(PC10 = Idents(Ang1.PC10),PC15 = Idents(Ang1.PC15))

# Doublet detection (scDblFinder)
# Convert Seurat object to SingleCellExperiment
sce.doublets <- as.SingleCellExperiment(Ang1.PC15)

# Inspect doublets in clusters and samples 
table(sce.doublets$orig.ident)
table(sce.doublets$seurat_clusters)

# Run scDblFinder
set.seed(123)
sce.doublets <- scDblFinder(sce.doublets,samples = "orig.ident",clusters = "seurat_clusters")

Ang1.PC15$doublet_class <- colData(sce.doublets)$scDblFinder.class
Ang1.PC15$doublet_score <- colData(sce.doublets)$scDblFinder.score

# Canonical markers from the reference paper
canonical_markers <- c(
  "Ms4a1",   # B cells
  "Ncr1",    # NK cells
  "Cd3e",    # T cells
  "S100a9",  # Granulocytes
  "Cd209a",  # DC-like cells
  "Fcgr1",   # Macrophages
  "Msln",    # Epicardial cells
  "Npr3",    # Endocardial cells
  "Pecam1",  # Endothelial cells
  "Lyve1",   # Lymphatic endothelial cells
  "Acta2",   # Smooth muscle cells
  "Pdgfrb",  # Pericytes
  "Pdgfra",  # Fibroblasts
  "Ttn",     # Cardiomyocytes
  "Plp1"     # Schwann cells
)

# BEFORE doublet removal: canonical marker DotPlot
canonical_present_before <- intersect(
  canonical_markers,
  rownames(Ang1.PC15)
)

DotPlot(
  Ang1.PC15,
  features = canonical_present_before,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(
    title = "Canonical cell-type markers before doublet removal",
    x = "Marker gene",
    y = "Cluster"
  )

# Visualize predicted doublets
DimPlot(Ang1.PC15,reduction = "umap",group.by = "doublet_class")

FeaturePlot(Ang1.PC15,features = "doublet_score")
table(Ang1.PC15$doublet_class)

prop.table(table(Ang1.PC15$seurat_clusters,Ang1.PC15$doublet_class),margin = 1)

# Before doublet removal: UMAP coloured by treatment
DimPlot(Ang1.PC15,reduction = "umap",group.by = "treatment")

# Before doublet removal
table(Ang1.PC15$treatment)

# Marker genes before doublet removal 
Ang1.PC15.markers <- FindAllMarkers(
  Ang1.PC15,
  only.pos = TRUE
)

head(Ang1.PC15.markers, n = 5)

# Select top 10 markers per cluster by avg_log2FC
top_markers_before <- Ang1.PC15.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

top_markers_before %>% filter(cluster == 0)
top_markers_before %>% filter(cluster == 1)
top_markers_before %>% filter(cluster == 2)
top_markers_before %>% filter(cluster == 3)
top_markers_before %>% filter(cluster == 4)
top_markers_before %>% filter(cluster == 5)
top_markers_before %>% filter(cluster == 6)
top_markers_before %>% filter(cluster == 7)
top_markers_before %>% filter(cluster == 8)
top_markers_before %>% filter(cluster == 9)
top_markers_before %>% filter(cluster == 10)
top_markers_before %>% filter(cluster == 11)
top_markers_before %>% filter(cluster == 12)

# Remove predicted doublets
# Remove cells predicted to be doublets and keep only singlets.
Ang1.PC15.filtered <- subset(Ang1.PC15,subset = doublet_class == "singlet")

# Insepct cells after doublet removal 
cat("After doublet removal:",
    "Genes =", nrow(Ang1.PC15.filtered),
    "Cells =", ncol(Ang1.PC15.filtered), "\n")

# Re-run normalization and dimensional reduction
Ang1.PC15.filtered <- NormalizeData(Ang1.PC15.filtered,normalization.method = "LogNormalize",scale.factor = 10000)
# Feature selection
Ang1.PC15.filtered <- FindVariableFeatures(Ang1.PC15.filtered,selection.method = "vst",nfeatures = 2000)
# Scaling
all.genes <- rownames(Ang1.PC15.filtered)
Ang1.PC15.filtered <- ScaleData(Ang1.PC15.filtered,features = all.genes)

# PCA again after doublet removal
Ang1.PC15.filtered <- RunPCA(Ang1.PC15.filtered,features = VariableFeatures(Ang1.PC15.filtered))

# Check PCA (elbow)
ElbowPlot(Ang1.PC15.filtered,ndims = 30)

# Clustering after doublet removal
Ang1.PC15.filtered <- FindNeighbors(Ang1.PC15.filtered,dims = 1:15)

# UMAP after doublet removal 
set.seed(123)
Ang1.PC15.filtered <- RunUMAP(Ang1.PC15.filtered,dims = 1:15)

# Test different clustering resolutions
# Compare several graph-clustering resolutions to assess how strongly the dataset partitions into finer cell populations.
for (res in c(0.2, 0.3,0.4, 0.5, 0.6, 0.7)) {Ang1.PC15.filtered <- FindClusters(Ang1.PC15.filtered,resolution = res)
print(DimPlot(Ang1.PC15.filtered,reduction = "umap",group.by = paste0("RNA_snn_res.", res),label = TRUE) + ggtitle(paste("Resolution", res)))}

# set the resolution 
Idents(Ang1.PC15.filtered) <- "RNA_snn_res.0.3"
Ang1.PC15.filtered$final_cluster <-Idents(Ang1.PC15.filtered) 


# AFTER doublet removal: canonical marker DotPlot
canonical_present_after <- intersect(
  canonical_markers,
  rownames(Ang1.PC15.filtered)
)

DotPlot(
  Ang1.PC15.filtered,
  features = canonical_present_after,
  group.by = "final_cluster"
) +
  RotatedAxis() +
  labs(
    title = "Canonical cell-type markers after doublet removal",
    x = "Marker gene",
    y = "Final cluster"
  )

# AFTER doublet removal: published marker FeaturePlots
# Published marker FeaturePlots
FeaturePlot(Ang1.PC15.filtered, features = c("Upk3b", "Gm12840", "Msln", "Gpc3", "Fmod"), ncol = 3)          # Epicardial
FeaturePlot(Ang1.PC15.filtered, features = c("Cytl1", "Vwf", "Ptgs1", "Cgnl1", "Plvap"), ncol = 3)           # Endocardial
FeaturePlot(Ang1.PC15.filtered, features = c("Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d"), ncol = 3)             # Cardiomyocytes
FeaturePlot(Ang1.PC15.filtered, features = c("Col1a1", "Crispld2", "Ogn", "Islr", "Pdgfra"), ncol = 3)       # Fibroblasts
FeaturePlot(Ang1.PC15.filtered, features = c("Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1"), ncol = 3)        # Endothelial
FeaturePlot(Ang1.PC15.filtered, features = c("Mmrn1", "Ccl21a", "Lyve1", "Cldn5", "Flt4"), ncol = 3)         # Lymphatic endothelial
FeaturePlot(Ang1.PC15.filtered, features = c("Tagln", "Myh11", "Olfr558", "Lmod1", "Nrip2"), ncol = 3)       # Smooth muscle
FeaturePlot(Ang1.PC15.filtered, features = c("Vtn", "Colec11", "Steap4", "Kcnj8", "Abcc9"), ncol = 3)        # Pericytes
FeaturePlot(Ang1.PC15.filtered, features = c("Kcna1", "Scn7a", "Plp1", "Gfra3", "Gpr37l1"), ncol = 3)        # Schwann
FeaturePlot(Ang1.PC15.filtered, features = c("Cd79a", "Ly6d", "H2-DMb2", "Cd79b", "Ms4a1"), ncol = 3)        # B cells
FeaturePlot(Ang1.PC15.filtered, features = c("H2-Q7", "Hcst", "Gimap3", "Il7r", "Skap1"), ncol = 3)          # T/NK cells
FeaturePlot(Ang1.PC15.filtered, features = c("Lmnb1", "Slpi", "Retnlg", "Clec4d", "Il1r2"), ncol = 3)        # Granulocytes
FeaturePlot(Ang1.PC15.filtered, features = c("Csf1r", "Adgre1", "Pld4", "Ms4a6c", "Mgl2"), ncol = 3)         # Macrophages

# After doublet removal: UMAP coloured by treatment
DimPlot(Ang1.PC15.filtered,reduction = "umap",group.by = "treatment")
table(Ang1.PC15.filtered$treatment)

DimPlot(Ang1.PC15.filtered,reduction = "umap",group.by = "final_cluster",label = TRUE)

# Compare treatment distribution across clusters
# Number of cells per cluster by treatment
table(Cluster = Ang1.PC15.filtered$final_cluster,Treatment = Ang1.PC15.filtered$treatment)

cluster_treatment_counts <- table(Cluster = Ang1.PC15.filtered$final_cluster,Treatment = Ang1.PC15.filtered$treatment)

prop.table(cluster_treatment_counts,margin = 2)

Ang1.PC15.filtered.markers <- FindAllMarkers(
  Ang1.PC15.filtered,
  only.pos = TRUE
)

# Select top 10 markers per cluster by avg_log2FC
top_markers_final <- Ang1.PC15.filtered.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

# View markers per cluster
top_markers_final %>% filter(cluster == 0)
top_markers_final %>% filter(cluster == 1)
top_markers_final %>% filter(cluster == 2)
top_markers_final %>% filter(cluster == 3)
top_markers_final %>% filter(cluster == 4)
top_markers_final %>% filter(cluster == 5)
top_markers_final %>% filter(cluster == 6)
top_markers_final %>% filter(cluster == 7)
top_markers_final %>% filter(cluster == 8)
top_markers_final %>% filter(cluster == 9)
top_markers_final %>% filter(cluster == 10)
top_markers_final %>% filter(cluster == 11)
top_markers_final %>% filter(cluster == 12)

# Heatmap of final markers
DoHeatmap(Ang1.PC15.filtered,features = top_markers_final$gene) + NoLegend()

top5_genes <- top_markers_final %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5
  ) %>%
  pull(gene) %>%
  unique()

avg_matrix <- AverageExpression(
  Ang1.PC15.filtered,
  features = top5_genes,
  group.by = "final_cluster",
  assays = "RNA"
)$RNA

# Rename columns (remove g prefix if present)
colnames(avg_matrix) <- gsub("^g", "", colnames(avg_matrix))

# Clean heatmap
pheatmap(avg_matrix,scale = "row",cluster_rows = FALSE,cluster_cols = FALSE, fontsize_row = 8, fontsize_col = 10)

# CREATE A SEPARATE OBJECT FOR RESOLUTION 0.7
Ang1.PC15.filtered.07 <- Ang1.PC15.filtered

Idents(Ang1.PC15.filtered.07) <- "RNA_snn_res.0.7"
Ang1.PC15.filtered.07$final_cluster_07 <- Idents(Ang1.PC15.filtered.07)

# UMAP AT RESOLUTION 0.7
DimPlot(
  Ang1.PC15.filtered.07,
  reduction = "umap",
  group.by = "final_cluster_07",
  label = TRUE
) +
  labs(title = "UMAP at resolution 0.7")

# Compare old 0.3 clusters with new 0.7 clusters
table(
  Resolution_0.3 = Ang1.PC15.filtered$final_cluster,
  Resolution_0.7 = Ang1.PC15.filtered.07$final_cluster_07
)
# CANONICAL MARKER DOTPLOT
canonical_present_07 <- intersect(
  canonical_markers,
  rownames(Ang1.PC15.filtered.07)
)

DotPlot(
  Ang1.PC15.filtered.07,
  features = canonical_present_07,
  group.by = "final_cluster_07"
) +
  RotatedAxis() +
  labs(
    title = "Canonical cell-type markers at resolution 0.7",
    x = "Marker gene",
    y = "Cluster"
  )

# FIND MARKER GENES FOR RESOLUTION 0.7 CLUSTERS
# Identify positively enriched genes for each resolution 0.7 cluster and rank markers by average log2 fold change.
Idents(Ang1.PC15.filtered.07) <- "final_cluster_07"

markers_07 <- FindAllMarkers(
  Ang1.PC15.filtered.07,
  only.pos = TRUE
)

# Select top 10 markers per cluster
top_markers_07 <- markers_07 %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

# Print top markers for every cluster
for (i in unique(top_markers_07$cluster)) {
  
  cat("\nCluster", i, "\n")
  
  print(
    top_markers_07 %>%
      filter(cluster == i)
  )
}

# MARKER HEATMAP
DoHeatmap(
  Ang1.PC15.filtered.07,
  features = unique(top_markers_07$gene),
  group.by = "final_cluster_07"
) +
  NoLegend()

# CLEAN AVERAGE-EXPRESSION HEATMAP
top5_genes_07 <- top_markers_07 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  pull(gene) %>%
  unique()

avg_matrix_07 <- AverageExpression(
  Ang1.PC15.filtered.07,
  features = top5_genes_07,
  group.by = "final_cluster_07",
  assays = "RNA"
)$RNA
colnames(avg_matrix_07) <- gsub(
  "^g",
  "",
  colnames(avg_matrix_07)
)
pheatmap(
  avg_matrix_07,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top cluster markers at resolution 0.7"
)

DotPlot(
  Ang1.PC15.filtered.07,
  features = c(
    # Epicardial
    "Upk3b", "Gm12840", "Msln", "Gpc3", "Fmod",
    
    # Endocardial
    "Cytl1", "Vwf", "Ptgs1", "Cgnl1", "Plvap",
    
    # Fibroblasts
    "Col1a1", "Crispld2", "Ogn", "Islr", "Pdgfra",
    
    # Cardiomyocytes
    "Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d",
    "Tnnc3", "Tnnt2", "Acta1", "Actc1",
    
    # Endothelial
    "Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1",
    "Pecam1", "Flt1", "Egfl7", "Tie1", "Epas", "Fabp4", 
    
    # Lymphatic endothelial
    "Mmrn1", "Ccl21a", "Lyve1", "Cldn5", "Flt4",
    
    # Smooth muscle
    "Tagln", "Myh11", "Olfr558", "Lmod1", "Nrip2",
    
    # Pericytes
    "Vtn", "Colec11", "Steap4", "Kcnj8", "Abcc9",
    
    # Schwann
    "Kcna1", "Scn7a", "Plp1", "Gfra3", "Gpr37l1",
    
    # B cells
    "Cd79a", "Ly6d", "H2-DMb2", "Cd79b", "Ms4a1",
    
    # T/NK cells
    "H2-Q7", "Hcst", "Gimap3", "Il7r", "Skap1",
    
    # Granulocytes
    "Lmnb1", "Slpi", "Retnlg", "Clec4d", "Il1r2",
    
    # Macrophages
    "Csf1r", "Adgre1", "Pld4", "Ms4a6c", "Mgl2"
  ),
  group.by = "final_cluster_07"
) +
  RotatedAxis()

DotPlot(
  Ang1.PC15.filtered,
  features = c(
    # Epicardial
    "Upk3b", "Gm12840", "Msln", "Gpc3", "Fmod",
    
    # Endocardial
    "Cytl1", "Vwf", "Ptgs1", "Cgnl1", "Plvap",
    
    # Cardiomyocytes
    "Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d",
    
    # Endothelial
    "Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1",
    
    # Fibroblasts
    "Col1a1", "Crispld2", "Ogn", "Islr", "Pdgfra",
    
    # Lymphatic endothelial
    "Mmrn1", "Ccl21a", "Lyve1", "Cldn5", "Flt4",
    
    # Smooth muscle
    "Tagln", "Myh11", "Olfr558", "Lmod1", "Nrip2",
    
    # Pericytes
    "Vtn", "Colec11", "Steap4", "Kcnj8", "Abcc9",
    
    # Schwann
    "Kcna1", "Scn7a", "Plp1", "Gfra3", "Gpr37l1",
    
    # B cells
    "Cd79a", "Ly6d", "H2-DMb2", "Cd79b", "Ms4a1",
    
    # T/NK cells
    "H2-Q7", "Hcst", "Gimap3", "Il7r", "Skap1",
    
    # Granulocytes
    "Lmnb1", "Slpi", "Retnlg", "Clec4d", "Il1r2",
    
    # Macrophages
    "Csf1r", "Adgre1", "Pld4", "Ms4a6c", "Mgl2"
  ),
  group.by = "final_cluster"
) +
  RotatedAxis()

focused_07 <- subset(
  Ang1.PC15.filtered.07,
  subset = final_cluster_07 %in% c("2", "5", "16")
)

markers_07 <- list(
  Cardiomyocyte = c("Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d", "Actn1", "Tnnt2", "Acta1", "Tnnc3"),
  Endothelial = c("Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1", "Pecam1", "Flt1", "Egfl7", "Tie1", "Epas", "Fabp4")
)

DotPlot(
  focused_07,
  features = markers_07,
  group.by = "final_cluster_07"
) +
  RotatedAxis() +
  labs(
    title = "Higher-resolution annotation of cluster 3",
    x = "Marker gene",
    y = "Resolution 0.7 cluster"
  )

# CELL-TYPE ANNOTATION
celltype_map_03 <- c(
  "0"  = "Endothelial",
  "1"  = "Endothelial",
  "2"  = "Fibroblast",
  "3"  = "Cardiomyocyte/endothelial",
  "4"  = "Macrophage",
  "5"  = "Smooth muscle",
  "6"  = "Pericyte",
  "7"  = "Endocardial",
  "8"  = "Fibroblast",
  "9"  = "Schwann",
  "10" = "B cell",
  "11" = "T/NK cell",
  "12" = "Endothelial"
)

Ang1.PC15.filtered$cell_type <- unname(
  celltype_map_03[Ang1.PC15.filtered$final_cluster]
)

# CELL-TYPE ANNOTATION UMAP
DimPlot(
  Ang1.PC15.filtered,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "Cell-type annotation",
    x = "UMAP 1",
    y = "UMAP 2"
  )

# RESOLUTION 0.7: INVESTIGATION OF CLUSTER 3
# Keep resolution 0.7 clusters derived from cluster 3
focused_07 <- subset(
  Ang1.PC15.filtered.07,
  subset = final_cluster_07 %in% c("2", "5", "16")
)

# Cardiomyocyte and endothelial marker sets
markers_07_focus <- list(
  Cardiomyocyte = c("Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d", "Actn2", "Tnnt2", "Acta1", "Actc1", "Tnnc3"),
  Endothelial = c("Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1", "Pecam1", "Flt1", "Egfl7", "Tie1", "Epas1", "Fabp4")
)

# Set cluster order
focused_07$final_cluster_07 <- factor(
  focused_07$final_cluster_07,
  levels = c("2", "5", "16")
)

# Cardiomyocyte and endothelial marker DotPlot
DotPlot(
  focused_07,
  features = markers_07_focus,
  group.by = "final_cluster_07",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "Cardiomyocyte and endothelial marker expression",
    x = "Marker gene",
    y = "Resolution 0.7 cluster"
  )

# DIFFERENTIAL EXPRESSION ANALYSIS
# Ang II vs Saline within resolution 0.3 clusters
# Create cluster-treatment identities
Ang1.PC15.filtered$cluster_treatment <- paste(
  Ang1.PC15.filtered$final_cluster,
  Ang1.PC15.filtered$treatment,
  sep = "_"
)

Idents(Ang1.PC15.filtered) <- "cluster_treatment"

# Run DE within each cluster
# Test Ang II versus saline expression separately within every resolution 0.3 cluster.
DE_results <- list()

for (i in levels(Ang1.PC15.filtered$final_cluster)) {
  
  DE_results[[paste0("Cluster_", i)]] <- FindMarkers(
    Ang1.PC15.filtered,
    ident.1 = paste0(i, "_Ang"),
    ident.2 = paste0(i, "_Sal"),
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

# Keep significant treatment-responsive genes
# Retain genes passing the adjusted-P-value and absolute log2-fold-change thresholds, then rank by effect size.
DE_results_filtered <- lapply(
  DE_results,
  function(x) {
    
    x %>%
      filter(
        p_val_adj < 0.05,
        abs(avg_log2FC) > 1
      ) %>%
      arrange(desc(abs(avg_log2FC)))
  }
)

# Number of significant DE genes per cluster
DE_gene_counts <- data.frame(
  Cluster = names(DE_results_filtered),
  Significant_DE_genes = sapply(
    DE_results_filtered,
    nrow
  )
)
print(DE_gene_counts)

# Display top 10 significant genes per cluster
for (i in names(DE_results_filtered)) {
  
  cat("\n", i, "\n")
  
  print(
    head(
      DE_results_filtered[[i]],
      10
    )
  )
}

# DIFFERENTIAL EXPRESSION ANALYSIS
# Resolution 0.7 clusters 2, 5 and 16
# Keep clusters derived from resolution 0.3 cluster 3
# Restrict the higher-resolution object to subclusters 2, 5, and 16 for focused treatment comparisons.
cluster3_07 <- subset(
  Ang1.PC15.filtered.07,
  subset = final_cluster_07 %in% c("2", "5", "16")
)

# Check cell numbers by cluster and treatment
table(
  Cluster_07 = cluster3_07$final_cluster_07,
  Treatment = cluster3_07$treatment
)

# Create cluster-treatment identities
cluster3_07$cluster_treatment_07 <- paste(
  cluster3_07$final_cluster_07,
  cluster3_07$treatment,
  sep = "_"
)

Idents(cluster3_07) <- "cluster_treatment_07"

# Run DE within each resolution 0.7 cluster
DE_cluster3_07 <- list()

for (i in c("2", "5", "16")) {
  
  DE_cluster3_07[[paste0("Cluster_", i)]] <- FindMarkers(
    cluster3_07,
    ident.1 = paste0(i, "_Ang"),
    ident.2 = paste0(i, "_Sal"),
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

# Keep significant treatment-responsive genes
DE_cluster3_07_filtered <- lapply(
  DE_cluster3_07,
  function(x) {
    
    x %>%
      filter(
        p_val_adj < 0.05,
        abs(avg_log2FC) > 1
      ) %>%
      arrange(desc(abs(avg_log2FC)))
  }
)

# Number of significant DE genes
DE_cluster3_07_counts <- data.frame(
  Cluster_07 = names(DE_cluster3_07_filtered),
  Significant_DE_genes = sapply(
    DE_cluster3_07_filtered,
    nrow
  )
)

print(DE_cluster3_07_counts)

# Display top 10 significant genes
for (i in names(DE_cluster3_07_filtered)) {
  
  cat("\n", i, "\n")
  
  print(
    head(
      DE_cluster3_07_filtered[[i]],
      10
    )
  )
}
