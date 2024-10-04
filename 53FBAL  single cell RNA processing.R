############# 53FBAL  single cell RNA processing

library(Seurat)
library("scDblFinder")
BiocManager::install("scDblFinder")
install.packages("NMF")
install.packages("iterators")
install.packages("pheatmap")

library("pheatmap")
iterators

# Define the paths to your CellRanger results
paths <- list(
  NT = ".../53FBAL_sc/NT/filtered_feature_bc_matrix",
  Tam = ".../Tam/filtered_feature_bc_matrix",
  Vand = ".../53FBAL_sc/Vand/outs/filtered_feature_bc_matrix"
)

# Create Seurat objects for each sample and add cell identifiers
seurat_objects <- lapply(names(paths), function(sample) {
  data <- Read10X(data.dir = paths[[sample]])
  seurat <- CreateSeuratObject(counts = data, project = sample)
  seurat <- RenameCells(seurat, new.names = paste(sample, Cells(seurat), sep = "_"))
  return(seurat)
})

# Merge the Seurat objects into one
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Save the merged Seurat object
saveRDS(merged_seurat, file = "merged_seurat.rds")

# Print the merged Seurat object to check the result
print(merged_seurat)

Assays(merged_seurat)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

merged_seurat

head(merged_seurat)
merged_seurat@assays


###########################

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

# List of layer names
layers <- c("counts.NT.SeuratProject", "counts.Tam.SeuratProject", "counts.Vand")
doublet_labels_list <- list()

for (layer in layers) {
  # Extract the specific layer
  layer_data <- GetAssayData(merged_seurat, assay = "RNA", slot = layer)
  
  # Create a new Seurat object with the extracted layer
  new_seurat <- CreateSeuratObject(counts = layer_data)
  
  # Convert to SingleCellExperiment and run scDblFinder
  new_seurat.sce <- as.SingleCellExperiment(new_seurat)
  new_seurat.sce <- scDblFinder(new_seurat.sce, clusters=FALSE)
  
  # Extract doublet labels
  doublet_labels <- new_seurat.sce$scDblFinder.class
  
  # Ensure doublet labels are named by cell barcodes
  names(doublet_labels) <- colnames(new_seurat.sce)
  
  # Store the doublet labels in a list
  doublet_labels_list[[layer]] <- doublet_labels
}

# Initialize a vector for doublet labels for all cells
all_doublet_labels <- rep("singlet", ncol(merged_seurat))
names(all_doublet_labels) <- colnames(merged_seurat)

# Combine doublet labels for all layers
for (layer in layers) {
  doublet_labels <- doublet_labels_list[[layer]]
  cell_indices <- which(names(all_doublet_labels) %in% names(doublet_labels))
  all_doublet_labels[cell_indices] <- doublet_labels[names(all_doublet_labels)[cell_indices]]
}

# Add doublet labels to the metadata of the original Seurat object
merged_seurat <- AddMetaData(merged_seurat, metadata = all_doublet_labels, col.name = "doublet_labels")


# Identify non-doublet cells
non_doublets <- merged_seurat$doublet_labels != "doublet"

# Subset the Seurat object to keep only non-doublet cells
merged_seurat_non_doublets <- subset(merged_seurat, cells = colnames(merged_seurat)[non_doublets])

# Verify the subset
table(merged_seurat_non_doublets$doublet_labels)

head(merged_seurat_non_doublets)


# Identify non-doublet cells
non_doublets <- merged_seurat$doublet_labels == "singlet"

# Subset the Seurat object to keep only non-doublet cells
merged_seurat_non_doublets <- subset(merged_seurat, cells = colnames(merged_seurat)[non_doublets])

# Verify the subset
table(merged_seurat_non_doublets$doublet_labels)


# Identify non-doublet cells
non_doublets <- merged_seurat$doublet_labels == 1

# Subset the Seurat object to keep only non-doublet cells
merged_seurat_non_doublets <- subset(merged_seurat, cells = colnames(merged_seurat)[non_doublets])

# Verify the subset
table(merged_seurat_non_doublets$doublet_labels)


merged_seurat_non_doublets_filtered <- subset(merged_seurat_non_doublets, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

merged_seurat_non_doublets_filtered <- SCTransform(merged_seurat_non_doublets_filtered, vars.to.regress = "percent.mt", verbose = FALSE)

DefaultAssay(merged_seurat_non_doublets_filtered)


### CEll annotation with SingleR

### SingleR cell type annotation


BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scRNAseq")


library(SingleR)
library(celldex)
library(Seurat)
library(scRNAseq)


library(scRNAseq)

fluidigm <- ReprocessedFluidigmData()
fluidigm
.libPaths("C:/library")


# Check for package updates: Ensure that all your packages, especially GenomeInfoDb, 
# GenomicRanges, SummarizedExperiment, and scRNAseq, are up to date.
# You can update packages using BiocManager::install()
# if you're using Bioconductor packages, or update.packages() for CRAN packages.

BiocManager::install("GenomeInfoDb")

library(GenomeInfoDb)

library(celldex)
ref <- HumanPrimaryCellAtlasData()
ref


library(SingleR)


test_assay <- GetAssayData(merged_seurat_non_doublets_filtered)
ref_assay <- GetAssayData(ref)

predictions <- SingleR(test = test_assay, ref = ref, labels = ref$label.main)

head(predictions)



# Extract the predicted labels from SingleR results
predicted_labels <- data.frame(cell = rownames(predictions), cell_type = predictions$pruned.labels)

# Ensure the rownames of the Seurat object match the cell identifiers in the predictions
rownames(predicted_labels) <- predicted_labels$cell

# Add the predicted cell types to the Seurat object metadata
merged_seurat_non_doublets_filtered <- AddMetaData(merged_seurat_non_doublets_filtered, metadata = predicted_labels$cell_type, col.name = "predicted_cell_type")

# Check to ensure the metadata has been added correctly
head(merged_seurat_non_doublets_filtered@meta.data)

# Save the modified Seurat object
saveRDS(merged_seurat_non_doublets_filtered, file = "merged_seurat_with_predictions.rds")
DefaultAssay(merged_seurat_non_doublets_filtered) <- "SCT"
# These are now standard steps in the Seurat workflow for visualization and clustering
merged_seurat_non_doublets_filtered <- RunPCA(merged_seurat_non_doublets_filtered, verbose = FALSE)
merged_seurat_non_doublets_filtered <- RunUMAP(merged_seurat_non_doublets_filtered, dims = 1:30, verbose = FALSE)

merged_seurat_non_doublets_filtered <- FindNeighbors(merged_seurat_non_doublets_filtered, dims = 1:30, verbose = FALSE)
merged_seurat_non_doublets_filtered <- FindClusters(merged_seurat_non_doublets_filtered, verbose = FALSE)
DimPlot(merged_seurat_non_doublets_filtered, label = TRUE)



