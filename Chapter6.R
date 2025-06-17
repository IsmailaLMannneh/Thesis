# Chapter 6 ------


#-----------------------------------------------------------------------------#
# PART 1: Geneeral workflow -----
#-----------------------------------------------------------------------------#

# i. Load required libraries ----
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggside)
library(scales)
library(pheatmap)
library(presto)
library(ggraph)
library(ExperimentHub)
library(harmony)
library(DESeq2)
library(cowplot)
library(scater)
library(dplyr)
library(tidyr)
library(stats)
library(ggplot2)
library(stringr)
library(multtest)
library(metap)
library(scDblFinder)
library(Azimuth)


# ii. Create data objects ------

# get data location 
dirs <- list.dirs(path = 'data/', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  print(paste("Processing directory:", x, "with name:", name))
  
  # read matrix files from the directories 
  cts <- ReadMtx(mtx = paste0('data/',x,'/matrix.mtx.gz'),
                 features = paste0('data/',x,'/features.tsv.gz'),
                 cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects 
  assign(name, CreateSeuratObject(counts = cts))
}

ls()
nameID = ls()[4:22]
nameID = nameID[-c(10)]


# merge datasets -
merged_seurat <- merge(NC08623_Converter, 
                       y = c(NC08744_Converter,NC08745_Converter,NC09158_Converter,
                             NC09181_Converter,NC09635_Resistor,NC09637_Resistor,
                             NC09639_Resistor,NC09641_Converter,NC09657_Converter,
                             NC09675_Converter,NC09679_Resistor, NC09685_Resistor,
                             NC09686_Resistor,NC09687_Resistor,NC09700_Resistor,
                             NC09702_Converter,NC09779_Resistor,NC09781_Converter),
                             add.cell.ids = nameID, project = 'NCProject')


# iii.  QC & filtering -----
merged_seurat$sample <- rownames(merged_seurat@meta.data)
# split sample column into sample id, group and cell barcode 
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('SampleID', 'Group', 'Barcode'), sep = '_')
# inspect the separation 
unique(merged_seurat@meta.data$Group) 
unique(merged_seurat@meta.data$SampleID)   
# Create mitochondria metadata and calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')


# View data by metrics 
# Visualize 
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
# individually by sample and group 
VlnPlot(merged_seurat, features = "nFeature_RNA", group.by = "SampleID", split.by = "Group", log = T)
VlnPlot(merged_seurat, features = "nCount_RNA", group.by = "SampleID", split.by = "Group", log = T)
VlnPlot(merged_seurat, features = "mitoPercent", group.by = "SampleID", split.by = "Group", log = T)
table(merged_seurat$Group) 

# Filtering poor quality cells: keep 
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 & 
                                   nFeature_RNA < 2500 & mitoPercent < 30)

# compare before and after
table(merged_seurat$Group) # Before cell and counts filter 
table(merged_seurat_filtered$Group) # After cell and counts filter 
table(merged_seurat_filtered$SampleID) 

# rename dataset as pbmc
pbmc = merged_seurat_filtered


# iv. Remove doublets  ------

pbmc_obj = merged_seurat_filtered
pbmc_obj = JoinLayers(pbmc_obj)
sce <- scDblFinder(GetAssayData(pbmc_obj, slot= "counts")) 

set.seed(1234)
results <- scDblFinder(sce, returnType = 'table') %>%
  as.data.frame() %>%
  filter(type == 'real')
head(results)

### 2030 cells (10.5%) doublets called 
keep = results %>%
  dplyr::filter(class == "singlet") %>%
  rownames()
pbmc_clean = pbmc_obj[, keep]
pbmc_clean

# v. Normalization using SCTrans
set.seed(12345)
pbmc <- SCTransform(object = pbmc_clean)
pbmc <- RunPCA(object = pbmc)
elbowplot = ElbowPlot(pbmc, 50)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE) # 30 looks best 
DimPlot(object = pbmc, reduction = 'umap')



# vi. Integrating SCTransformed data
# Batch correction}
pbmc <- IntegrateLayers(object = pbmc, method = HarmonyIntegration,
                         normalization.method = "SCT",
                         verbose = F)

pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc) # default res=0.8 gives 23 comms


# vii. run Azimuth cell anotations ----

# Load reference data
pbmcref<-readRDS("./Reference/ref.Rds")
pbmc = JoinLayers(pbmc)

# Load Query data
pbmc_lab <- RunAzimuth(pbmc, reference = "./Reference/") # Requires Internet


# vii. Exploring different resolutions using cluster tree ----
pbmc <- FindClusters(object = pbmc, resolution = 
                        c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
library(clustree)
scPlot <- clustree(pbmc, prefix="SCT_snn_res.")
DimPlot(pbmc, reduction = 'umap', group.by = c('Group',"SCT_snn_res.0.1", "SampleID"))
DimPlot(pbmc, reduction = 'umap', group.by = c("SCT_snn_res.0.1"))



#-----------------------------------------------------------------------------#
# PART 2: Cell proportion comparisons -----
#-----------------------------------------------------------------------------#

library(dplyr)
library(Seurat)

# Step 1: Extract the cell type labels from the assay
celltypes <- GetAssayData(pbmc, assay = "prediction.score.celltype.l1", slot = "data")
celltype_labels <- colnames(celltypes)[apply(celltypes, 2, which.max)]  

# Step 2: Add cell type labels to metadata
pbmc$CellType_L1 <- celltype_labels

# Step 3: Create dataframe with cell type, SampleID, Group
meta <- pbmc@meta.data %>%
  select(SampleID, Group, CellType_L1)

# Step 4: Count number of each cell type per sample
cell_counts <- meta %>%
  group_by(SampleID, Group, CellType_L1) %>%
  summarise(count = n(), .groups = "drop")

# Step 5: Total cells per sample
total_counts <- meta %>%
  group_by(SampleID) %>%
  summarise(total = n(), .groups = "drop")

# Step 6: Merge and compute proportions
cell_props <- left_join(cell_counts, total_counts, by = "SampleID") %>%
  mutate(proportion = count / total)

# Step 7: Run Wilcoxon tests per cell type
wilcox_results <- cell_props %>%
  group_by(CellType_L1) %>%
  summarise(
    p.value = wilcox.test(proportion ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

print(wilcox_results)



#-----------------------------------------------------------------------------#
# PART 3: DEG ANALYSES -----
#-----------------------------------------------------------------------------#

# i. CD4 T cells 
pbmc = PrepSCTFindMarkers(pbmc)
CD4T_RST.Markers <- FindMarkers(pbmc, ident.1 = "Resistor_CD4 T", ident.2 = "Converter_CD4 T", 
                                assay = "SCT", min.pct = 0, 
                                logfc.threshold = 0)
CD4T_RST.Markers <- tibble::rownames_to_column(CD4T_RST.Markers, "gene")
dim(CD4T_RST.Markers) 

# Filter results by p values and fold change for significant genes 
CD4T_RST_upregulated <- CD4T_RST.Markers[CD4T_RST.Markers$avg_log2FC >=0.25 & 
                                           CD4T_RST.Markers$p_val_adj < 0.05,"gene"] 
CD4T_RST_downregulated <- CD4T_RST.Markers[CD4T_RST.Markers$avg_log2FC <=-0.25 
                                           & CD4T_RST.Markers$p_val_adj < 0.05,"gene"]

DEG_CD4_combine = rbind(CD4T_RST_upregulated, CD4T_RST_downregulated)
write.csv(DEG_CD4_combine, "DEG_CD4_combine.csv")


# ii. B cells

B_RST.Markers <- FindMarkers(pbmc, ident.1 = "Resistor_B", ident.2 = "Converter_B", 
                             assay = "SCT", min.pct = 0, logfc.threshold = 0)
B_RST.Markers <- tibble::rownames_to_column(B_RST.Markers, "gene")

# Filter by p values and fold change 
B_RST_upregulated <- B_RST.Markers[B_RST.Markers$avg_log2FC >=0.25 & 
                                     B_RST.Markers$p_val_adj < 0.05,"gene"]
B_RST_downregulated <- B_RST.Markers[B_RST.Markers$avg_log2FC <=-0.25 
                                     & B_RST.Markers$p_val_adj < 0.05,"gene"]

DEG_BCells_combine = rbind(B_RST_upregulated, B_RST_downregulated)
write.csv(DEG_BCells_combine, "DEG_BCells_combine.csv")


# iii. CD8 T cells
CD8T_RST.Markers <- FindMarkers(pbmc, ident.1 = "Resistor_CD8 T", ident.2 = "Converter_CD8 T", 
                                assay = "SCT", min.pct = 0, logfc.threshold = 0)
CD8T_RST.Markers <- tibble::rownames_to_column(CD8T_RST.Markers, "gene")

# Filter by p values and fold change 
CD8T_RST_upregulated <- CD8T_RST.Markers[CD8T_RST.Markers$avg_log2FC >=0.25 & 
                                           CD8T_RST.Markers$p_val_adj < 0.05,"gene"]
CD8T_RST_downregulated <- CD8T_RST.Markers[CD8T_RST.Markers$avg_log2FC <=-0.25 
                                           & CD8T_RST.Markers$p_val_adj < 0.05,"gene"]

DEG_CD8Cells_combine = rbind(CD8T_RST_upregulated, CD8T_RST_downregulated)
write.csv(DEG_CD8Cells_combine, "DEG_CD8Cells_combine.csv")


# iv. Monocyte
Mono_RST.Marker <- FindMarkers(pbmc, ident.1 = "Resistor_Mono", ident.2 = "Converter_Mono", 
                               assay = "SCT", min.pct = 0, logfc.threshold = 0)
Mono_RST.Marker <- tibble::rownames_to_column(Mono_RST.Marker, "gene")
# Filter by p values and fold change 
Mono_RST_upregulated <- Mono_RST.Marker[Mono_RST.Marker$avg_log2FC >=0.25 & 
                                          Mono_RST.Marker$p_val_adj < 0.05,"gene"] 
Mono_RST_downregulated <- Mono_RST.Marker[Mono_RST.Marker$avg_log2FC <=-0.25 
                                          & Mono_RST.Marker$p_val_adj < 0.05,"gene"] 

DEG_Monocytes_combine = rbind(Mono_RST_upregulated, Mono_RST_downregulated)
write.csv(DEG_Monocytes_combine, "DEG_Monocytes_combine.csv")



# v. NK comparison
NK_RST.Marker <- FindMarkers(pbmc, ident.1 = "Resistor_NK", ident.2 = "Converter_NK", 
                             assay = "SCT", min.pct = 0, logfc.threshold = 0)
NK_RST.Marker <- tibble::rownames_to_column(NK_RST.Marker, "gene")
# Filter by p values and fold change 
NK_RST_upregulated <- NK_RST.Marker[NK_RST.Marker$avg_log2FC >=0.25 & 
                                      NK_RST.Marker$p_val_adj < 0.05,"gene"] 
NK_RST_downregulated <- NK_RST.Marker[NK_RST.Marker$avg_log2FC <=-0.25 
                                      & NK_RST.Marker$p_val_adj < 0.05,"gene"] 

DEG_NKcells_combine = rbind(NK_RST_upregulated, NK_RST_downregulated)
write.csv(DEG_NKcells_combine, "DEG_NKcells_combine.csv")



# vi. other T
OtherT_RST.Marker <- FindMarkers(pbmc, ident.1 = "Resistor_other T", ident.2 = "Converter_other T", 
                                 assay = "SCT", min.pct = 0, logfc.threshold = 0)
OtherT_RST.Marker <- tibble::rownames_to_column(OtherT_RST.Marker, "gene")
# Filter by p values and fold change 
OtherT_RST_upregulated <- OtherT_RST.Marker[OtherT_RST.Marker$avg_log2FC >=0.25 & 
                                              OtherT_RST.Marker$p_val_adj < 0.05,"gene"] 
OtherT_RST_downregulated <- OtherT_RST.Marker[OtherT_RST.Marker$avg_log2FC <=-0.25 
                                              & OtherT_RST.Marker$p_val_adj < 0.05,"gene"] 

DEG_OtherT_cells_combine = rbind(OtherT_RST_upregulated, OtherT_RST_downregulated)
write.csv(DEG_OtherT_cells_combine, "DEG_OtherT_cells_combine.csv")



# vii. DC comparison ######### No difference
DC_RST.Marker <- FindMarkers(pbmc, ident.1 = "Resistor_DC", ident.2 = "Converter_DC", 
                             assay = "SCT", min.pct = 0, logfc.threshold = 0)
DC_RST.Marker <- tibble::rownames_to_column(DC_RST.Marker, "gene")
# Filter by p values and fold change 
DC_RST_upregulated <- DC_RST.Marker[DC_RST.Marker$avg_log2FC >=0.25 & 
                                      DC_RST.Marker$p_val_adj < 0.05,"gene"] 
DC_RST_downregulated <- DC_RST.Marker[DC_RST.Marker$avg_log2FC <=-0.25 
                                      & DC_RST.Marker$p_val_adj < 0.05,"gene"] 

DEG_DC_cells_combine = rbind(DC_RST_upregulated, DC_RST_downregulated)
write.csv(DEG_DC_cells_combine, "DEG_DC_cells_combine.csv")


# viii. other comparison
OtherCells_RST.Marker <- FindMarkers(pbmc, ident.1 = "Resistor_other", ident.2 = "Converter_other", 
                                     assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25)
OtherCells_RST.Marker <- tibble::rownames_to_column(OtherCells_RST.Marker, "gene")
# Filter by p values and fold change 
OtherCells_RST_upregulated <- OtherCells_RST.Marker[OtherCells_RST.Marker$avg_log2FC >=0.25 & 
                                                      OtherCells_RST.Marker$p_val_adj < 0.05,"gene"] 
OtherCells_RST_downregulated <- OtherCells_RST.Marker[OtherCells_RST.Marker$avg_log2FC <=-0.25 
                                                      & OtherCells_RST.Marker$p_val_adj < 0.05,"gene"]


DEG_other_cells_combine = rbind(OtherCells_RST_upregulated, OtherCells_RST_downregulated)
write.csv(DEG_other_cells_combine, "DEG_other_cells_combine.csv")


# Combine all DEGs data frames
combined_DE_genes <- bind_rows(
  DEG_CD4_combine,
  DEG_BCells_combine,
  DEG_CD8Cells_combine,
  DEG_Monocytes_combine,
  DEG_NKcells_combine,
  DEG_OtherT_cells_combine,
  DEG_DC_cells_combine,
  DEG_other_cells_combine,
  .id = "Cell_Type")

# View first few rows of the combined data
head(combined_DE_genes)
write.csv(combined_DE_genes, "combined_DE_genes_new.csv")




#-----------------------------------------------------------------------------#
# PART 4: Plots, GSEA, Cell chat and gene modules scores ANALYSES -----
#-----------------------------------------------------------------------------#


####### PLOTS ------

############## Figure Panel 1: QC and Clustering -------

# A. QC plot -----
#-------------------------------------------------------------#
# use pic in slides 


# B. Integrated data by Sample ID ----
#-------------------------------------------------------------#
pbmc
# reduction= pbmc@reductions[["integrated_dr"]]
## Umap
pdf("Panel 1B_UMAP for Integrated data_march2025.pdf",width = 8, height = 6)
set.seed(112)
DimPlot(pbmc, reduction = "umap", label = F, group.by = 'SampleID')
dev.off()



# C. Seurat Clusters ----
#-------------------------------------------------------------#

####Reoder levels
levels(pbmc)
# my_levels = c( "0" , "1","2" , "3" , "4" , "5" , "6" , "7" , "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
# combined.dataset@meta.data$seurat_clusters = factor(x=combined.dataset@meta.data$seurat_clusters, levels = my_levels)

## Umap
pdf("Panel 1_UMAP_for_Suerat_clusters_march2025.pdf",width = 8, height = 6)
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = 'seurat_clusters')
dev.off()


# D. labelled Seurat Clusters -----
#-------------------------------------------------------------#

Idents(pbmc) = "predicted.celltype.l1"

## Umap
pdf("Panel 1_UMAP_for_Li_Azimuth_clusters_march2025.pdf",width = 8, height = 6)
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = 'predicted.celltype.l1')
dev.off()


## Umap split 
pdf("Panel 1_UMAP_for_Li_Azimuth_clusters_Plit by group march2025.pdf",width = 10, height = 6)
DimPlot(pbmc, reduction = "umap", label = F, group.by = 'predicted.celltype.l1', split.by = "Group")
dev.off()


# E. Heat map of cluster markers -----
#-------------------------------------------------------------#

Idents(pbmc) = "seurat_clusters"

####Finding Marker genes
Combined.dataset_markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.dataset_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(x=Combined.dataset_markers, file = "Marker genes for all 21 clusters_march2025.csv", quote = FALSE)

### A heatmap of the top 10 marker genes
top10 = Combined.dataset_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

Heatmap10 = DoHeatmap(pbmc, features = top10$gene, size = 8, angle =0)+ scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(axis.text.y = element_text(size = 15))+
  theme(axis.text = element_text(face="bold"))
Heatmap10 = ggsave("Panel 1. Heatmap for top 10 marker genes_in_n21 clusters_march2025.jpg",width = 25, height = 25, units = "in", dpi = 300)


# Heat map of top5 only 
### A heatmap of the top 10 marker genes
top5 = Combined.dataset_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

Heatmap5 = DoHeatmap(pbmc, features = top5$gene, size = 8, angle =0)+ scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(axis.text.y = element_text(size = 15))+
  theme(axis.text = element_text(face="bold"))

Heatmap5 = ggsave("Panel 1_Heatmap for top 5 marker genes_in_n21 clusters_march2025.jpg",width = 25, height = 25, units = "in", dpi = 300)



# F. Dot plot of cluster specific markers -------
#-------------------------------------------------------------#

# Plot top2 genes  

top2 = Combined.dataset_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2


DPlot1 = DotPlot(pbmc, features = unique(top2$gene[!top2$gene %in% c("MT-CO3", "MT-ND4")]), 
                 group.by = "seurat_clusters", cols = c("blue", "red")) + 
  RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)) + 
  xlab(NULL) + ylab(NULL) + ggtitle("Lineage-specific Markers") + 
  theme(plot.title = element_text(hjust = 0.5))

pdf("Panel_1_Dot_plot_top_markers.pdf",width = 8, height = 6)
DPlot1
dev.off()

# Lineage specific TFs  

lineage_TFs <- c("TBX21", "GATA3", "FOXP3", "RORC",  # T cells  
                 "PAX5", "BCL6", "PRDM1", "IRF4",  # B cells  
                 "SPI1", "IRF8", "NFKB1", "CEBPA",  # Monocytes  
                 "EOMES", "NFIL3", "ZEB2",  # NK cells  
                 "BATF3", "ZBTB46")  # Dendritic cells  


# Check which genes from lineage_TFs are present
genes_found <- lineage_TFs %in% rownames(pbmc@assays$SCT@counts)
names(genes_found) <- lineage_TFs
genes_found

# Dotplot 2 
DPlot2 = DotPlot(pbmc, features = lineage_TFs, 
                 group.by = "seurat_clusters", cols = c("blue", "red")) + 
  RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Lineage-specific TFs") + theme(plot.title = element_text(hjust = 0.5))
DPlot2 # Not great at all. 

### Plot the two
pdf("Panel 1 Marker and TF factor plots.pdf", height = 10, width = 10)
DPlot1/DPlot2
dev.off()

###### 
markers_to_plot = c(cell_markers <- c(
  "CD19", "CD20", "PAX5",    # B cells
  "CD3", "CD4", "CD8",    # T cells
  "CD56", "NKG2D", "NKp46",    # NK cells
  "CD11c", "HLA-DR", "CLEC10A",  # Dendritic cells
  "CD14", "CD68", "CD163")    # Macrophages/Monocytes
)

# Dotplot 2 
DPlot3 = DotPlot(pbmc, features = markers_to_plot, 
                 group.by = "seurat_clusters", cols = c("blue", "red")) + 
  RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Lineage-specific TFs") + theme(plot.title = element_text(hjust = 0.5))
DPlot3




##----------------------------------------------------------##
# Figure Panel 2: Cell composition and DEG------
##----------------------------------------------------------##

# D. Table of DEGS

table(combined_DE_genes2$CellID)

# Create summary table
cell_summary <- as.data.frame(table(combined_DE_genes2$CellID))

# Check output
print(cell_summary)

# Save as CSV
write.csv(cell_summary, "Figure_T_Table_CellID_Summary.csv", row.names = FALSE)



# E. Volcano plot ------

library(ggplot2)
library(ggrepel)

# Define significance threshold
logFC_threshold <- 0.25
pval_threshold <- 0.05

# Volcano plot
ggplot(combined_DE_genes2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC > logFC_threshold)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  facet_wrap(~ CellID) +  # Facet for each cell type
  scale_color_manual(values = c("blue", "red")) +  # Blue for downregulated, red for upregulated
  labs(x = "log2 Fold Change", y = "-log10 Adjusted P-value", title = "Volcano Plots of DEGs in Different Cell Types") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(data = combined_DE_genes2 %>% filter(p_val_adj < pval_threshold & abs(avg_log2FC) > logFC_threshold),
                  aes(label = gene), size = 3)

ggsave("Figure A. Volcano plots in one.pdf", height = 8, width = 10)


#### Amended Volcano plot 
library(dplyr)
library(ggplot2)
library(ggrepel)

# Define thresholds
logFC_threshold <- 0.25
pval_threshold <- 0.05

# Get top 5 up and top 5 downregulated genes per cell type
top_genes <- combined_DE_genes2 %>%
  filter(p_val_adj < pval_threshold) %>%
  group_by(CellID) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 4) %>%  # Top 5 upregulated
  bind_rows(
    combined_DE_genes2 %>%
      filter(p_val_adj < pval_threshold) %>%
      group_by(CellID) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 4))  # Top 5 downregulated

write.csv(top_genes, "top_genes.csv")

# Amend Force label to make top_gene2 in excel. and load 
library(readr)
top_genes2 <- read_csv("top_genes2.csv")
# overwrite 
top_genes = top_genes2


top_genes$CellID
top_genes = top_genes[1:45, ]

# Volcano plot
ggplot(combined_DE_genes2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC > logFC_threshold)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray") +
  facet_wrap(~ CellID) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "log2 Fold Change", y = "-log10 Adjusted P-value", title = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 13)
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 100
  )

ggsave("Figure 6.4. Volcano plots in one.pdf", height = 8, width = 10)



# F. Violin like plot ------

library(ggplot2)
library(ggrepel)

# Define thresholds for labeling top DE genes
logFC_threshold <- 1  
pval_threshold <- 0.05

# Select top 6 DE genes per CellID based on absolute Log2FC
top_DE_genes <- combined_DE_genes2 %>%
  filter(p_val_adj < pval_threshold) %>%
  group_by(CellID) %>%
  slice_max(order_by = abs(avg_log2FC), n = 6)  # Top 5 per CellID

# Create violin plot
ggplot(combined_DE_genes2, aes(x = CellID, y = avg_log2FC, fill = CellID)) +
  geom_jitter(aes(color = avg_log2FC > 0), width = 0.08, alpha = 0.3, size = 2) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add baseline at 0
  geom_text_repel(data = top_DE_genes, aes(label = gene), 
                  size = 3, max.overlaps = 25, force = 2) + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +  # Nice color scheme
  labs(x = "Cell Type", y = "Log2 Fold Change", 
       title = "Violin Plot of Log2 Fold Change per Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure B. Violin like plot DEGS.pdf", height = 8, width = 10)




########  presention for figure panle 2 cell proportions: -----

## cell number L1 

clusters <- data.frame(table(pbmc$predicted.celltype.l1, pbmc$Group)) #gives you a data frame of raw counts of the meta.data you are interested in each cluster - can change the resolution as needed
colnames(clusters) <- c("Cluster", "Group", "Count") #rename the columns

clusters2 <- data.frame(table(pbmc$Group)) # Totals
colnames(clusters2) <- c("Group", "Total") #rename the columns

clusters3 <- left_join(clusters, clusters2, by = "Group") #combine the data frames
clusters3$Frequency <- clusters3$Count/clusters3$Total #create proportion category (can multiple by 100 if want percent). 

#DimPlot(pbmc, group.by = "predicted.celltype.l1", label = T, , label.size = 4) + NoLegend() |
clusters3 %>% 
  ggplot(aes(x = Cluster, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("darkred", "darkblue"))

ggsave("Figure_panel 2_Umap and Cell numbers L1.pdf", height = 6, width = 10)


## cell numbers L2

clusters <- data.frame(table(pbmc$predicted.celltype.l2, pbmc$Group)) #gives you a data frame of raw counts of the meta.data you are interested in each cluster - can change the resolution as needed
colnames(clusters) <- c("Cluster", "Group", "Count") #rename the columns

clusters2 <- data.frame(table(pbmc$Group)) # Totals
colnames(clusters2) <- c("Group", "Total") #rename the columns

clusters3 <- left_join(clusters, clusters2, by = "Group") #combine the data frames
clusters3$Frequency <- clusters3$Count/clusters3$Total #create proportion category (can multiple by 100 if want percent). 

#DimPlot(pbmc, group.by = "predicted.celltype.l2", label = T, label.size = 4) + NoLegend() |


#### Exclude cells found in only one group 
clusters3 <- clusters3 %>% filter(!Cluster %in% c("CD8 Proliferating", "HSPC"))

clusters3 %>% 
  ggplot(aes(x = Cluster, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  #scale_fill_manual(values = c("darkred", "darkblue"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure_6.2.pdf", height = 4, width = 6)

ggsave("Figure_panel 2_Umap and Cell numbers L1.pdf", height = 6, width = 14)



# Note code: Voiline plot for individual markers -----
Vln1 = VlnPlot(obj = pbmc, features = 'CD3D', group.by = "predicted.celltype.l1") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 
Vln2 = VlnPlot(obj = pbmc, features = 'KLRC1', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Vln3 = VlnPlot(obj = pbmc, features = 'SCGB3A1', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())

## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1, Vln2, Vln3, ncol = 2)

# Combined_Plots = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12)

# Violins with p values ----

Violn1 = VlnPlot(obj = pbmc, features = 'CD3D', group.by = "Group", cols = c("red", "blue"), y.max = 7) + #scale_x_discrete(limits=c("WT", "HDAC1cKO")) 
  theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + 
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("Resistor", "Converter", method = "wilcox.test")))
Violn1 



##----------------------------------------------------------##
# Figure Panel 4: Cellchat ------
##----------------------------------------------------------##

### Part 1: Set up individula cell chat objects ------
#----------------------------------------------------------##

options(future.globals.maxSize= 8*1024^3)

# Cell chat part1 (setting up the chat objects)

# Load required libraries
library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(NMF)


# Set up CellChat object for Resistors
pbmc_resistor <- subset(pbmc, subset = Group == "Resistor") # ,assay = "RNA"
cellchat_resistor <- createCellChat(object = pbmc_resistor, group.by = "predicted.celltype.l1")

# Set up CellChat object for Converters
pbmc_converter <- subset(pbmc, subset = Group == "Converter")
cellchat_converter <- createCellChat(object = pbmc_converter, group.by = "predicted.celltype.l1")


#remove R objects to save memory
rm(list = setdiff(ls(), c("pbmc_resistor", "pbmc_converter",
                          "cellchat_converter", "cellchat_resistor")))

# Load CellChatDB human database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat_resistor@DB <- CellChatDB
cellchat_converter@DB <- CellChatDB


# Preprocessing (Identify overexpressed genes and interactions)
cellchat_resistor <- subsetData(cellchat_resistor) # This step is necessary even if using the whole database
cellchat_resistor <- identifyOverExpressedGenes(cellchat_resistor)
cellchat_resistor <- identifyOverExpressedInteractions(cellchat_resistor)
#
cellchat_converter <- subsetData(cellchat_converter)
cellchat_converter <- identifyOverExpressedGenes(cellchat_converter)
cellchat_converter <- identifyOverExpressedInteractions(cellchat_converter)

# Compute communication probability
ptm = Sys.time()
cellchat_resistor <- computeCommunProb(cellchat_resistor) #default is Trimean, other Trim (more intx)
cellchat_resistor <- filterCommunication(cellchat_resistor) #  min.cells = 10, can filter
cellchat_converter <- computeCommunProb(cellchat_converter)
cellchat_converter <- filterCommunication(cellchat_converter)

# Extract the inferred cellular communication network as a data frame
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest.

df.net_R<- subsetCommunication(cellchat_resistor)
df.net_C <- subsetCommunication(cellchat_converter)


# Compute aggregated communication networks
cellchat_resistor <- computeCommunProbPathway(cellchat_resistor)
cellchat_resistor <- aggregateNet(cellchat_resistor)
cellchat_converter <- computeCommunProbPathway(cellchat_converter)
cellchat_converter <- aggregateNet(cellchat_converter)

# Compare the number of interactions and interaction strength
#par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_resistor@net$count,
                        vertex.weight = 0.01,
                        weight.scale = T,
                        label.edge= F,
                        title.name = "Number of interactions")
#gg1 <- netVisual_circle(cellchat_resistor@net$weight,
#          vertex.weight = groupSize,
#          title.name = "Resistor", title.name = "Number of interactions")
# vertex.weight=groupSize (instead of 0.01), weight.scale = T, label.edge= F,

gg2 <- netVisual_circle(cellchat_converter@net$count,
                        vertex.weight = 0.01,
                        weight.scale = T,
                        label.edge= F,
                        title.name = "Number of interactions")
#gg2 <- netVisual_circle(cellchat_converter@net$count,
#                        vertex.weight = groupSize,
#                        title.name = "Converter", title.name = "Interaction strength")

library(patchwork)
gg1 | gg2

# Visualize individual cell communications
# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.

# For resistors
png("Resistors_individual_cell_communications_panel.png",
    width = 10, height = 10, res = 1200, units = "in")

mat <- cellchat_resistor@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = 0.01, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


# For converters
png("Converter_individual_cell_communications_panel.png",
    width = 10, height = 10, res = 1200, units = "in")

mat <- cellchat_converter@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = 0.01, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Note: interpretations: Thicker edge line indicates a stronger signal.
# In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram.

# Visualization options : heirachy and chord
netVisual_chord_cell
netVisual_individual

# Visualization of cell-cell communication at different levels
netVisual_aggregate

# All the signaling pathways showing significant communications can be accessed by
cellchat_resistor@netP$pathways
cellchat_converter@netP$pathways

# # Save results
# save.image("CellChatObjects.RData")

#### Part 2: Continuation to Merge Object ---------

# Load objects from Part 1
load("CellChatObjects.RData")

#remove R objects to save memory 
rm(list = setdiff(ls(), c("pbmc_resistor", "pbmc_converter", 
                          "cellchat_converter", "cellchat_resistor")))

# Merge CellChat objects
object.list <- list(Converter = cellchat_converter, Resistor = cellchat_resistor)  
cellchatmerged <- mergeCellChat(object.list, add.names = names(object.list))
View(cellchatmerged)


# A. compare interactions

gg1 <- compareInteractions(cellchatmerged, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchatmerged, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("Figure Panel 3A Diffeenrce in total interactions resisters vs converters.pdf",
       height = 6, width = 8)


## B. Compare the number of interactions and interaction strength among different cell populations

p1 = netVisual_diffInteraction(cellchatmerged, weight.scale = T)
# red colored edges represent increased (blue, decreased) signaling in the second dataset (converters) compared to the first
p1
ggsave("Figure Panel 3B1 viz of in total interactions resisters vs converters count.pdf",
       height = 6, width = 8)
p2 = netVisual_diffInteraction(cellchatmerged, weight.scale = T, measure = "weight")
p2
ggsave("Figure Panel 3B2 viz of in total interactions resisters vs converters weight.pdf",
       height = 6, width = 8)



# C. Heatmap showing differential number of interactions or interaction strength

heatmap1 <- netVisual_heatmap(cellchatmerged)
heatmap1
ggsave("Figure Panel 3C1 Heatmap of total interactions.pdf",
       height = 8, width = 8)

heatmap2 <- netVisual_heatmap(cellchatmerged, measure = "weight")
heatmap2
ggsave("Figure Panel 3C1 Heatmap of interactions strength.pdf",
       height = 8, width = 8)

# View side by side 
heatmap1 + heatmap2


### (C) Circle plot 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
ggsave("Figure Panel 3B1 Alternative side by side comparisons.pdf",
       height = 6, width = 8)



# We then can show the number of interactions or interaction strength between any two cell types in each dataset. 
# {r, fig.width=9,fig.height = 4, fig.wide = TRUE, fig.align = "center"}
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Identify cell populations with significant changes in sending or receiving signals between different datasets by following option A, and the signaling changes of specific cell populations by following option B. 

### (A) Identify cell populations with significant changes in sending or receiving signals
# {r, fig.width=9,fig.height = 4, fig.wide = TRUE, fig.align = "center"}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

ggsave("Figure pane 3D_ incoming and ougoing signals.pdf", height = 6, width = 8)

#> Results: From the scatter plot, we can see that Monocytes and DCs emerge as one of the major source and targets in resistors 
#> compared to conveerters.


### (B) Identify the signaling changes of specific cell populations

# Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. 
# 
# cellchatmerged Mono and DC
cellsig1 <- netAnalysis_signalingChanges_scatter(cellchatmerged, idents.use = "Mono", signaling.exclude = "MIF")
cellsig2 <- netAnalysis_signalingChanges_scatter(cellchatmerged, idents.use = "DC", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(cellsig1,cellsig2))

ggsave("Figure panel 4_monocyte differences.pdf", height = 6, width = 8)



# #> Part II: Identify altered signaling 
cellchatmerged <- computeNetSimilarityPairwise(cellchatmerged, type = "functional")
cellchatmerged <- netEmbedding(cellchatmerged, type = "functional")
cellchatmerged <- netClustering(cellchatmerged, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchatmerged, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchatmerged, type = "functional", nCol = 2)

### Identify signaling groups based on structure similarity
cellchatmerged <- computeNetSimilarityPairwise(cellchatmerged, type = "structural")
cellchatmerged <- netEmbedding(cellchatmerged, type = "structural")
cellchatmerged <- netClustering(cellchatmerged, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchatmerged, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchatmerged, type = "structural", nCol = 2)


### (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
# CellChat can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network). 
# 
# This bar chart can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. When setting `do.stat = TRUE`, a paired Wilcoxon test is performed to determine whether there is a significant difference of the signaling information flow between two conditions. The top signaling pathways colored red are enriched in NL skin, and these colored greens are enriched in the LS skin.

# {r, fig.width=9,fig.height = 4, fig.wide = TRUE, fig.align = "center"}
rankplot1 <- rankNet(cellchatmerged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
rankplot2 <- rankNet(cellchatmerged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

rankplot1 + rankplot2

ggsave("Figure panel 3_ Signalling pathways.pdf", height = 6, width = 10)

### (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ggsave("Figure Panle 3E Outgoing signal side by side.pdf", height = 6, width = 8)


# Incoming
hti1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
hti2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(hti1 + hti2, ht_gap = unit(0.5, "cm"))
ggsave("Figure Panle 3E incoming signal side by side.pdf", height = 6, width = 8)

# All signals 
htAll1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
htAll2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(htAll1 + htAll2, ht_gap = unit(0.5, "cm"))
ggsave("Figure Panle 3E All signal side by side.pdf", height = 6, width = 8)



##----------------------------------------------------------##
# Panel 3: GSEA ------
##----------------------------------------------------------##

# # Load required libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(pheatmap)

# Load your dataset
df <- read.csv("combined_DE_genes_with_cellID.csv")

# Ensure column names match your dataset
head(df)

# Load MSigDB Hallmark pathways
pathways <- msigdbr(species = "Homo sapiens", category = "H") 
pathways_list <- split(pathways$gene_symbol, pathways$gs_name)

# Create output folder
dir.create("GSEA_results", showWarnings = FALSE)

# Get unique cell types
cell_types <- unique(df$CellID)

# Function to perform GSEA for each cell type
run_gsea <- function(cell_type, df, pathways_list) {
  
  # Filter data for the specific cell type
  cell_data <- df %>% filter(CellID == cell_type) %>% select(gene, avg_log2FC)
  
  # Rank genes by log fold change
  ranks <- setNames(cell_data$avg_log2FC, cell_data$gene)
  
  # Run fgseaMultilevel
  fgseaRes <- fgseaMultilevel(pathways = pathways_list, stats = ranks)
  
  # Sort by significance
  fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
  
  # Convert list columns (e.g., leadingEdge) to strings
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  
  # Save results
  write.csv(fgseaRes, paste0("GSEA_results/", gsub(" ", "_", cell_type), "_GSEA_results.csv"), row.names = FALSE)
  
  # Generate plot for the top enriched pathway
  if (nrow(fgseaRes) > 0) {
    topPathway <- fgseaRes$pathway[1]
    p <- plotEnrichment(pathways_list[[topPathway]], ranks) + 
      ggtitle(paste("GSEA for", cell_type, "\nTop Pathway:", topPathway))
    
    # Save plot
    ggsave(paste0("GSEA_results/", gsub(" ", "_", cell_type), "_GSEA_plot.png"), p, width = 6, height = 4)
  }
}

# Run GSEA for all cell types
for (cell_type in cell_types) {
  run_gsea(cell_type, df, pathways_list)
}



########## visualization:

# Set directory where GSEA results are stored
result_dir <- "GSEA_results"

# Get list of all GSEA result files
files <- list.files(result_dir, pattern = "_GSEA_results.csv", full.names = TRUE)

# Initialize an empty list to store data
gsea_list <- list()

# Read each file and extract NES values
for (file in files) {
  cell_type <- gsub("_GSEA_results.csv", "", basename(file))  # Extract cell type name
  gsea_data <- read.csv(file) %>%
    select(pathway, NES) %>%
    rename(!!cell_type := NES)  # Rename NES column based on cell type
  
  gsea_list[[cell_type]] <- gsea_data
}

# Merge data from all cell types into one dataframe
gsea_df <- Reduce(function(x, y) full_join(x, y, by = "pathway"), gsea_list)

# Convert to matrix for heatmap (remove pathway column)
rownames(gsea_df) <- gsea_df$pathway
gsea_matrix <- as.matrix(gsea_df[, -1])

# Replace NA values with 0 (if some pathways are missing in certain cell types)
gsea_matrix[is.na(gsea_matrix)] <- 0


# Dot plot :  -------

gsea_long <- gsea_df %>%
  pivot_longer(cols = -pathway, names_to = "Cell_Type", values_to = "NES") %>%
  mutate(Significance = ifelse(NES > 1, "Up", "Down"))  

# Remove Pathway and underscore 
gsea_long$pathway = gsub("HALLMARK_", "", gsea_long$pathway)
gsea_long$Cell_Type = gsub("_", "", gsea_long$Cell_Type)


# Dot Plot
pdf("Figure Panel 3 GSEA using Dot plot.pdf", height = 8, width = 8)
ggplot(gsea_long, aes(x = Cell_Type, y = pathway, size = abs(NES), color = Significance)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(title = "Pathway Enrichment Across Cell Types",
       x = "Cell Type",
       y = "Pathway",
       size = "NES (Absolute)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



##----------------------------------------------------------##
# Panel 5: Signature scores ------
##----------------------------------------------------------##

# set up memory 
options(future.globals.maxSize= 8*1024^3)

# load Pbmc from environ 
pbmc = readRDS("pbmc.RDS")

# Cytotoxicity module : Add the module scores to your Seurat object 
signature_genes_cytox <- c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM') 
pbmc <- AddModuleScore(pbmc, features = list(signature_genes_cytox), name = "CYTX_Signature")

# add score 

## ADD SIGNATURE SCORE TO THE METADATA COLUMNS
metadata <-pbmc@meta.data %>%
  mutate(Cytotocicity = 
           case_when(MySignature1 > 0.25 ~
                       "Cytotoxic",
                     MySignature1 <= 0.25 ~
                       "Non_cytotoxic"))
pbmc <- AddMetaData(pbmc, metadata)



# Inflammation Modules: Add to seurat object 
signature_genes_inflamm <- c("TNF", "NFKB1", "IL1B", "NFKB2", "CXCL10", "SOCS3","IKBKB", "MAPK8", "IL6", "CXCL8", "TRAF6", "ISG15", "MX1", "OAS1", "IFNG", "IRF1") # custom

# Add the module scores to your Seurat object using AddModuleScore
pbmc <- AddModuleScore(pbmc, features = list(signature_genes_inflamm), name = "full_inflammation_Signature")


## ADD  SCORE TO METADATA 
metadata <-pbmc@meta.data %>%
  mutate(full_inflammatory_module = 
           case_when(full_inflammation_Signature1 > 0.08 ~
                       "Inflammatry_cells",
                     full_inflammation_Signature1 <= 0.08 ~
                       "Non_Inflammatry_cells"))
pbmc <- AddMetaData(pbmc, metadata)



# Cells subsets  -----
#T and NK cells subset
TCellsubset = subset(pbmc, CellType %in% c("CD4 T", "CD8 T", "NK", "other T"))

# B cells 
Bcellsubset = subset(pbmc, CellType %in% c("B"))

# Monocyte subset
Monosubset = subset(pbmc, CellType %in% c("Mono"))

# DC subset 
seuratObjD = subset(pbmc, CellType %in% c("DC"))

# Innate cell subset 
MonocyteDCsubset = subset(pbmc, CellType %in% c("Mono", "DC"))


# Figure panel 6 :---- 

# 1. cytoxicity dimplot -----

### PLOT OF SIGNATURE BY GROUP 
pdf("Fig6.A2.dimplot.pdf", width = 5, height = 4)
DimPlot(pbmc, group.by = "Cytotocicity", split.by = "Group") #ILM
dev.off()

# 2. By group bar ------
pdf("Fig6.A2b.pdf", width = 4, height = 4)
df %>%
  ggplot(aes(x = Group, y = Proportion, fill = Cytotocicity)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c("Cytotoxic" = "#F8766D", "Non_cytotoxic" = "#00BFC4"))
dev.off()


# Rename vars 
df$Group = gsub("Resistor", "Resister", df$Group)
df$Cytotoxicity = df$Cytotocicity
df$Cytotoxicity = gsub("Non_cytotoxic", "Non cytotoxic", df$Cytotoxicity)

# Flip 
my_levels = c("Non cytotoxic", "Cytotoxic")
df$Cytotoxicity <- factor(x = df$Cytotoxicity, levels = my_levels)

# By group  
pdf("Fig6.A2b.dimplot.pdf", width = 4, height = 3)

df %>%
  ggplot(aes(x = Group, y = Proportion, fill = Cytotoxicity)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Cytotoxic" = "#F8766D", "Non cytotoxic" = "#00BFC4")) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
dev.off()


# 3a Distribution across resisters -----
TCellResistors = subset(TCellsubset, Group=="Resistor")
unique(TCellResistors$Group)


# flip plot 

### PROPORTION OF SIGNATURE POS AND NEG CELLS 
df_tcell <- data.frame(table(TCellResistors$SampleID, TCellResistors$Cytotocicity))
colnames(df_tcell) = c("SampleID", "Cytotocicity", "Count")

df_tcell1 <- data.frame(table(TCellResistors$SampleID))
colnames(df_tcell1) = c("SampleID", "Total")

df_tcell <- left_join(df_tcell, df_tcell1, by = "SampleID")
df_tcell$Proportion <- df_tcell$Count/df_tcell$Total

df_tcell %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Cytotocicity)) +
  geom_bar(stat = "identity", position = "fill")

# Rename 
names(df_tcell)
df_tcell$Cytotoxicity = df_tcell$Cytotocicity
df_tcell$Cytotoxicity = gsub("_", " ", df_tcell$Cytotoxicity)

# Recode sample IDs 
library(dplyr)

df_tcell$SampleID <- recode(df_tcell$SampleID,
                            "NC09635" = "R1",
                            "NC09637" = "R2",
                            "NC09639" = "R3",
                            "NC09679" = "R4",
                            "NC09685" = "R5",
                            "NC09686" = "R6",
                            "NC09687" = "R7",
                            "NC09700" = "R8",
                            "NC09779" = "R9")

# flip plot 
my_levels = c("Non cytotoxic", "Cytotoxic")
df_tcell$Cytotoxicity <- factor(x = df_tcell$Cytotoxicity, levels = my_levels)


pdf("Fig6.1B2.pdf", width = 6, height = 4)
df_tcell %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Cytotoxicity)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Cytotoxic" = "#F8766D", "Non cytotoxic" = "#00BFC4")) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
dev.off()


# 3b Distribution across Converters -----

TCellConverters = subset(TCellsubset, Group=="Converter")
unique(TCellConverters$Group)


# flip plot 

### PROPORTION OF SIGNATURE POS AND NEG CELLS 
ddf_tcell <- data.frame(table(TCellConverters$SampleID, TCellConverters$Cytotocicity))
colnames(ddf_tcell) = c("SampleID", "Cytotocicity", "Count")

ddf_tcell1 <- data.frame(table(TCellConverters$SampleID))
colnames(ddf_tcell1) = c("SampleID", "Total")

ddf_tcell <- left_join(ddf_tcell, ddf_tcell1, by = "SampleID")
ddf_tcell$Proportion <- ddf_tcell$Count/ddf_tcell$Total

ddf_tcell %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Cytotocicity)) +
  geom_bar(stat = "identity", position = "fill")

# Rename 
names(ddf_tcell)
ddf_tcell$Cytotoxicity = ddf_tcell$Cytotocicity
ddf_tcell$Cytotoxicity = gsub("_", " ", ddf_tcell$Cytotoxicity)

# Recode sample IDs 
library(dplyr)


ddf_tcell$SampleID <- recode(ddf_tcell$SampleID,
                             "NC08623" = "C1",
                             "NC08744" = "C2",
                             "NC08745" = "C3",
                             "NC09158" = "C4",
                             "NC09181" = "C5",
                             "NC09641" = "C6",
                             "NC09675" = "C7",
                             "NC09702" = "C8",
                             "NC09781" = "C9")

# flip plot 
my_levels = c("Non cytotoxic", "Cytotoxic")
ddf_tcell$Cytotoxicity <- factor(x = ddf_tcell$Cytotoxicity, levels = my_levels)


pdf("Fig6.1F.pdf", width = 6, height = 4)
ddf_tcell %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Cytotoxicity)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Cytotoxic" = "#F8766D", "Non cytotoxic" = "#00BFC4")) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
dev.off()


# 4. Cytotoxicity statistical comparisons -----
library(ggstatsplot)
# R object is TCellsubset
unique(TCellsubset$Cytotocicity)
TCell_cytoxcells = subset(TCellsubset,Cytotocicity=="Cytotoxic")

# rename group 
TCell_cytoxcells$Group = gsub("Resistor", "Resister", TCell_cytoxcells$Group)

pdf("Fig6.1D.pdf", height = 4, width = 4)
ggbetweenstats(
  data  = TCell_cytoxcells@meta.data[sample(1:nrow(TCell_cytoxcells@meta.data), 1605),],
  x    = "Group", 
  y    = "CYTX_Signature1",
  title = "Cytotoxicity_module response",
  type  = "nonparametric",
  point.args = list(alpha = 0))
dev.off()

# wilcox 
wilcox.test(CYTX_Signature1 ~ Group, data = TCell_cytoxcells@meta.data)



## inflammation Signature plots -----

# 1. Dimplot ------

pdf("Fig.6.1a1.pdf", width = 6, height = 4)
DimPlot(pbmc_Inf, group.by = "full_inflammatory_module", split.by = "Group")
dev.off()


# 2 group ------

pdf("Fig.6.1a2.pdf", width = 4, height = 4)
dtfl %>%
  ggplot(aes(x = Group, y = Proportion, fill = full_inflammatory_module)) +
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = c("Inflammatry_cells" = "#F8766D", "Non_Inflammatry_cells" = "#00BFC4"))
dev.off()

# rename and re-plot 
dtfl$Group = gsub("Resistor", "Resister", dtfl$Group)
dtfl$Inflammation = dtfl$full_inflammatory_module
dtfl$Inflammation = gsub("Non_Inflammatry_cells", "Non Inflammatory", dtfl$Inflammation)

dtfl$Inflammation = gsub("Inflammatry_cells", "Inflammatory", dtfl$Inflammation)


# Flip 
my_levels = c("Non Inflammatory", "Inflammatory")
dtfl$Inflammation <- factor(x = dtfl$Inflammation, levels = my_levels)

# By group  
pdf("Fig.6.1a2.pdf", width = 4, height = 3)

dtfl %>%
  ggplot(aes(x = Group, y = Proportion, fill = Inflammation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Inflammatory" = "#F8766D", "Non Inflammatory" = "#00BFC4")) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

dev.off()


# 3. By individual  bar plot -----

pbmc_R = subset(pbmc_Inf, Group=="Resistor")
unique(pbmc_R$Group)


#### By samples ---

# subset Resisiters 
pbmc_R = subset(pbmc_Inf, Group=="Resistor")
pbmc_C = subset(pbmc_Inf, Group=="Converter")

# flip plot 

### 3A. Distribution among Resisters ----
df_pbmc_Inf <- data.frame(table(pbmc_R$SampleID, pbmc_R$full_inflammatory_module))
colnames(df_pbmc_Inf) = c("SampleID", "full_inflammatory_module", "Count")

df_pbmc_Inf1 <- data.frame(table(pbmc_R$SampleID))
colnames(df_pbmc_Inf1) = c("SampleID", "Total")

df_pbmc_Inf <- left_join(df_pbmc_Inf, df_pbmc_Inf1, by = "SampleID")
df_pbmc_Inf$Proportion <- df_pbmc_Inf$Count/df_pbmc_Inf$Total

df_pbmc_Inf %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = full_inflammatory_module)) +
  geom_bar(stat = "identity", position = "fill")

# flip plot 
my_levels = c("Non_Inflammatry_cells", "Inflammatry_cells")
df_pbmc_Inf$full_inflammatory_module <- factor(x = df_pbmc_Inf$full_inflammatory_module, levels = my_levels)

# Distribution among converters 
df_pbmc_Inf$SampleID <- recode(df_tcell$SampleID,
                               "NC09635" = "R1",
                               "NC09637" = "R2",
                               "NC09639" = "R3",
                               "NC09679" = "R4",
                               "NC09685" = "R5",
                               "NC09686" = "R6",
                               "NC09687" = "R7",
                               "NC09700" = "R8",
                               "NC09779" = "R9")

# renames 
names(df_pbmc_Inf)
names(df_pbmc_Inf)[2] = "Inflammation"
df_pbmc_Inf$Inflammation = gsub("_", " ", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub(" cells", " ", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub("Inflammatry ", "Inflammatory", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub("Non Inflammatry ", "Non Inflammatory", df_pbmc_Inf$Inflammation)


# Flip levels 
my_levels = c("Non Inflammatory", "Inflammatory")
df_pbmc_Inf$Inflammation <- factor(x = df_pbmc_Inf$Inflammation, levels = my_levels)


pdf("Fig.6.1a3.pdf", width = 6, height = 4)
df_pbmc_Inf %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Inflammation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Inflammatory" = "#F8766D", "Non Inflammatory" = "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
dev.off()


### 
### 3B. Distribution among Converters ----
df_pbmc_Inf <- data.frame(table(pbmc_C$SampleID, pbmc_C$full_inflammatory_module))
colnames(df_pbmc_Inf) = c("SampleID", "full_inflammatory_module", "Count")

df_pbmc_Inf1 <- data.frame(table(pbmc_C$SampleID))
colnames(df_pbmc_Inf1) = c("SampleID", "Total")

df_pbmc_Inf <- left_join(df_pbmc_Inf, df_pbmc_Inf1, by = "SampleID")
df_pbmc_Inf$Proportion <- df_pbmc_Inf$Count/df_pbmc_Inf$Total

df_pbmc_Inf %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = full_inflammatory_module)) +
  geom_bar(stat = "identity", position = "fill")

# flip plot 
my_levels = c("Non_Inflammatry_cells", "Inflammatry_cells")
df_pbmc_Inf$full_inflammatory_module <- factor(x = df_pbmc_Inf$full_inflammatory_module, levels = my_levels)

# Recode sample IDs
df_pbmc_Inf$SampleID <- recode(df_pbmc_Inf$SampleID,
                             "NC08623" = "C1",
                             "NC08744" = "C2",
                             "NC08745" = "C3",
                             "NC09158" = "C4",
                             "NC09181" = "C5",
                             "NC09641" = "C6",
                             "NC09675" = "C7",
                             "NC09702" = "C8",
                             "NC09781" = "C9")

# renames 
names(df_pbmc_Inf)
names(df_pbmc_Inf)[2] = "Inflammation"
df_pbmc_Inf$Inflammation = gsub("_", " ", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub(" cells", " ", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub("Inflammatry ", "Inflammatory", df_pbmc_Inf$Inflammation)
df_pbmc_Inf$Inflammation = gsub("Non Inflammatry ", "Non Inflammatory", df_pbmc_Inf$Inflammation)


# Flip levels 
my_levels = c("Non Inflammatory", "Inflammatory")
df_pbmc_Inf$Inflammation <- factor(x = df_pbmc_Inf$Inflammation, levels = my_levels)


pdf("Fig.6.1b.pdf", width = 6, height = 4)
df_pbmc_Inf %>%
  ggplot(aes(x = SampleID, y = Proportion, fill = Inflammation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Inflammatory" = "#F8766D", "Non Inflammatory" = "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
dev.off()


# 4. Inflammation stats comparisons -----
library(ggstatsplot)
pbmc_Inf

# subset inflammatory cells  
pbmc_Inflamms = subset(pbmc_Inf,full_inflammatory_module=="Inflammatry_cells")

# rename group 
pbmc_Inflamms$Group = gsub("Resistor", "Resister", pbmc_Inflamms$Group)


# Compare between groups 

pdf("Fig.6.1a4.pdf", height = 4, width = 5)

ggbetweenstats(
  data  = pbmc_Inflamms@meta.data[sample(1:nrow(pbmc_Inflamms@meta.data), 708),],
  x    = "Group",
  y    = "full_inflammation_Signature1",
  title = "Comparison of Inflammatory_response",
  type  = "nonparametric",
  point.args = list(alpha = 0))

dev.off()

# wilcox 
wilcox.test(full_inflammation_Signature1 ~ Group, data = pbmc_Inflamms@meta.data)



# save step 
save.image("Figure6panelsFinal.RData")
