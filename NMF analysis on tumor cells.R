##################### NMF analysis

genes_vand_signature <- read.csv("~/genes_vand_signature.csv", sep="")
genes <- genes_vand_signature$Gene # Assuming the column with gene names is named 'GeneName'
Tumor_cells_53FBAL_NT_Vand_subset_vand_genes <- subset(Tumor_cells_53FBAL_NT_Vand, features = genes)

head(Tumor_cells_53FBAL_NT_Vand_subset_vand_genes)

Assays(Tumor_cells_53FBAL_NT_Vand_subset_vand_genes)

DefaultAssay(Tumor_cells_53FBAL_NT_Vand_subset_vand_genes)<- "SCT"

# Extract the expression matrix
expression_matrix <- GetAssayData(Tumor_cells_53FBAL_NT_Vand_subset_vand_genes, assay = "SCT", slot = "data")


# Convert the expression matrix to a format suitable for NMF
expression_matrix <- as.matrix(expression_matrix)

head(expression_matrix)

# Remove rows with only null or NA values
expression_matrix <- expression_matrix[rowSums(is.na(expression_matrix)) < ncol(expression_matrix) &
                                         rowSums(expression_matrix == 0, na.rm = TRUE) < ncol(expression_matrix), ]

rank <- 8 # Adjust this based on your requirements

# Run NMF
nmf_result <- nmf(expression_matrix, rank, .options = 't')

library(NMF)

basis_matrix <- basis(nmf_result)
coefficient_matrix <- coef(nmf_result)



# Extract the cell metadata from the Seurat object
cell_metadata <- Tumor_cells_53FBAL_NT_Vand_subset_vand_genes@meta.data

# Extract the prefixes and barcodes
cell_barcodes_with_prefix <- rownames(cell_metadata)
cell_barcodes_without_prefix <- sub("^[^_]*_", "", cell_barcodes_with_prefix)
prefixes <- sub("_[^_]*$", "", cell_barcodes_with_prefix)
prefix_mapping <- setNames(prefixes, cell_barcodes_without_prefix)

# Convert the coefficient matrix to a data frame for easier manipulation
coefficient_df <- as.data.frame(t(coefficient_matrix))

rownames(coefficient_df)


# Check the new row names
head(coefficient_df)
tail(coefficient_df)
# Optionally, save the coefficient data frame with the updated barcodes
write.csv(coefficient_df, "coefficient_matrix_NMF_vand_NT.csv", row.names = TRUE)

library(Seurat)

# Add a new column with the sample identity (orig.ident) from the Seurat metadata
coefficient_df$sample_identity <- cell_metadata[rownames(coefficient_df), "orig.ident"]

# Check the first few rows to ensure the sample identity has been added correctly
head(coefficient_df)


# Filter the data by sample identity
NT_df <- coefficient_df %>% filter(sample_identity == "NT") %>% select(-sample_identity)

Vand_df <- coefficient_df %>% filter(sample_identity == "Vand") %>% select(-sample_identity)



# Define custom color limits
color_limits <- c(-0.02,0.05)

# Plot heatmaps with ordered NMF factors
pheatmap(as.matrix(NT_df), cluster_rows = TRUE, cluster_cols = FALSE, 
         scale = "none", display_numbers = FALSE, 
         main = "Heatmap for NT", 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = seq(color_limits[1], color_limits[2], length.out = 101))

pheatmap(as.matrix(Vand_df), cluster_rows = TRUE, cluster_cols = FALSE, 
         scale = "none", display_numbers = FALSE, 
         main = "Heatmap for Vand", 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = seq(color_limits[1], color_limits[2], length.out = 101))

# Remove the sample_identity column for calculations
data_matrix <- coefficient_df %>% select(-sample_identity)

# Calculate the mean for every NMF factor for each sample identity
mean_factors <- coefficient_df %>%
  group_by(sample_identity) %>%
  summarise(across(everything(), mean))

# Print the mean factors
print(mean_factors)

# Assign each cell to the NMF factor with the maximum value
max_factor <- data_matrix %>%
  mutate(max_factor = apply(., 1, which.max))

# Add sample_identity back to max_factor
max_factor <- cbind(max_factor, sample_identity = coefficient_df$sample_identity)

# Calculate the percentage of cells in each factor for each sample identity
percentage_factors <- max_factor %>%
  group_by(sample_identity, max_factor) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Print the percentage of cells in each factor
print(percentage_factors)





# Ensure the NMF factors (columns) are ordered from 1 to 8
ordered_columns <- paste0("V", 1:8)
coefficient_df <- coefficient_df[, c(ordered_columns, "sample_identity")]

# Assign each cell to the NMF factor with the maximum value
data_matrix <- coefficient_df %>% select(-sample_identity)
max_factor <- data_matrix %>%
  mutate(max_factor = apply(., 1, which.max))

# Add sample_identity back to max_factor
max_factor <- cbind(max_factor, sample_identity = coefficient_df$sample_identity)

# Create a list to store mean scores for each NMF factor
mean_scores_list <- list()

# Generate a dataframe for each NMF factor and calculate the mean for the corresponding NMF factor
for (i in 1:8) {
  df_filtered <- max_factor %>% filter(max_factor == i)
  mean_scores <- df_filtered %>%
    group_by(sample_identity) %>%
    summarise(mean_score = mean(get(paste0("V", i)), na.rm = TRUE)) %>%
    mutate(factor = i) %>%
    ungroup()
  mean_scores_list[[i]] <- mean_scores
}

# Combine the results into a single dataframe for better readability
combined_mean_scores <- bind_rows(mean_scores_list)

# Print the combined mean scores
print(combined_mean_scores)

head(max_factor)



# Prepare the metadata
metadata <- max_factor %>%
  select(max_factor, sample_identity)


# Ensure the rownames of metadata match the cell names in the Seurat object
rownames(metadata) <- rownames(max_factor)

# Add metadata to the Seurat object
Tumor_cells_53FBAL_NT_Vand <- AddMetaData(Tumor_cells_53FBAL_NT_Vand, metadata = metadata)

# Verify the metadata
head(Tumor_cells_53FBAL_NT_Vand@meta.data)

Idents(Tumor_cells_53FBAL_NT_Vand)<- "max_factor"

DefaultAssay(Tumor_cells_53FBAL_NT_Vand) <- "SCT"

Tumor_cells_53FBAL_NT_Vand <-PrepSCTFindMarkers(Tumor_cells_53FBAL_NT_Vand)

# Perform DE analysis using FindMarkers
de_results_NMF1_vs_rest <- FindMarkers(Tumor_cells_53FBAL_NT_Vand, ident.1 = "1", ident.2 = c("2", "3", "4","5","6", "7", "8"), group.by = "max_factor")

write.csv(de_results_NMF1_vs_rest, "de_results_NMF1_vs_rest.csv")


# Perform DE analysis using FindMarkers
de_results_NMF6_vs_rest <- FindMarkers(Tumor_cells_53FBAL_NT_Vand, ident.1 = "6", ident.2 = c("2", "3", "4","5","1", "7", "8"), group.by = "max_factor")

write.csv(de_results_NMF6_vs_rest, "de_results_NMF6_vs_rest.csv")


############## Scoring gene signatures


signatures <- list(
  RET_signature_Rasha_up = c("CHGA", "MUC2", "HMP19", "MUC6", "GPR26", "ANK1", "UNC13A", "CNR1", "PXDN", "CPLX2", "SLIT1", "FGF10", "SNAP25", "KCNK3", "MYT1", "SEZ6L", "WDR17", "NXPH1", "ASCL1", "RIMS2", "RUNDC3A", "GABRB3", "PVALB", "CHRNA2", "CHRNB2", "CEACAM6", "PTPRN2", "KCNJ3", "DNAJC6", "INSM1", "CEL", "LRP2", "SYP", "NRXN3", "NMNAT2", "GLDN", "FAM5B", "CEACAM5", "NTRK2", "CACNB2", "CD177", "LDB3"),
  RET_signature_Rasha_down = c("LEP", "C4orf7", "LBP", "FAM5C", "STAC2", "PLA2G2A", "MUC16", "PNMT", "ACTA1", "MIA", "DCD", "COL9A3", "SYT8", "TCN1", "S100A9", "CALML3", "DHRS2", "DKK1", "KRT81", "RBM20", "KLK10", "DES", "KLK7", "S100A7", "CALML5", "S100A8", "OBP2B", "KRT16", "CBLN2", "KLK11", "MUCL1", "SCGB3A1", "CRISP3", "UGT2B4", "CYP4F8", "MYBPC1", "UGT2B11", "TSIX", "KRT6A", "KRT13")
)


library(UCell)

Tumor_cells_53FBAL_NT_Vand <- AddModuleScore_UCell(Tumor_cells_53FBAL_NT_Vand, 
                                                   features=signatures, name=NULL)


# Create violin plots for the three signature scores stored in the metadata
VlnPlot(Tumor_cells_53FBAL_NT_Vand, 
        features = "RET_signature_Rasha_down", 
        group.by = "max_factor", 
        pt.size = 0)



# Extract the metadata including the signature scores and orig.ident
signature_data <- FetchData(Tumor_cells_53FBAL_NT_Vand, vars = c("RET_signature_Rasha_down",
                                                         "max_factor"))

head(signature_data)


# Subset data for factor 1 and factor 6
factor1_data <- signature_data[signature_data$max_factor == 1, "RET_signature_Rasha_down"]
factor6_data <- signature_data[signature_data$max_factor == 6, "RET_signature_Rasha_down"]

# Calculate mean and median for factor 1
mean_factor1 <- mean(factor1_data)
median_factor1 <- median(factor1_data)

# Calculate mean and median for factor 6
mean_factor6 <- mean(factor6_data)
median_factor6 <- median(factor6_data)

# Perform Wilcoxon test between factor 1 and factor 6
wilcox_test <- wilcox.test(factor1_data, factor6_data)

# Display results
mean_factor1
median_factor1
mean_factor6
median_factor6
wilcox_test



########


# Extract the metadata including the signature scores and orig.ident
signature_data <- FetchData(Tumor_cells_53FBAL_NT_Vand, vars = c("RET_signature_Rasha_up",
                                                                 "max_factor"))

head(signature_data)

RET_signature_Rasha_up = wilcox.test(ERBB2_GOBP ~ orig.ident, data = signature_data)

# Subset data for factor 1 and factor 6
factor1_data <- signature_data[signature_data$max_factor == 1, "RET_signature_Rasha_up"]
factor6_data <- signature_data[signature_data$max_factor == 6, "RET_signature_Rasha_up"]

# Calculate mean and median for factor 1
mean_factor1 <- mean(factor1_data)
median_factor1 <- median(factor1_data)

# Calculate mean and median for factor 6
mean_factor6 <- mean(factor6_data)
median_factor6 <- median(factor6_data)

# Perform Wilcoxon test between factor 1 and factor 6
wilcox_test <- wilcox.test(factor1_data, factor6_data)

# Display results
mean_factor1
median_factor1
mean_factor6
median_factor6
wilcox_test


############ T47D down signature
\
"MT-ATP6", "FMN1", "NDFIP1", "SRP9", "IL6ST", "MRFAP1", "NDUFC2", "TFRC", "FKBP4", "H3F3A", "THBS1", "RPS10", "PPM1K", "ATP5MC3", "MARCKSL1", "MYBL1", "SFPQ", "OLFM1", "HNRNPH3", "PSENEN", "PSMA2", "PCLO", "MINOS1", "B4GALT1", "KRT19", "TIMMDC1", "ARPC5", "PMVK", "ARL3", "MRPL13", "SDF2", "NDUFA8", "CDV3", "MRPL22", "RUNX1", "IGFBP4", "TSFM", "SDSL", "TPBG", "MT-CO3", "GABARAPL2", "GNG11", "LINC00992", "KRT15", "ERH", "TTC14", "AHSA1", "SERPINA6", "REEP5", "LDHA", "SGK3", "MT-CYB", "SLTM", "CA12", "DNAJA1", "TRAPPC4", "GREB1", "HNRNPM", "NDUFB8", "FKBP5", "EIF4A3", "CCT3", "SRSF7", "PFKFB3", "GOLM1", "MAGED2", "PRDX4", "TRPS1", "TMBIM4", "MT-ND2", "HNRNPA2B1", "PGR", "CA8", "XBP1", "SNRPB", "PRLR", "RBBP7", "VDAC3", "CCT5", "PDCD4", "HSP90AA1", "NMU", "DEGS2", "BNIP3", "UBB", "SLC40A1", "JPT1", "TUBA1A", "SPA17", "CALM3", "ENO1", "TUBB4B", "GATA3", "CALM2", "STC2", "H2AFV", "IGFBP5", "ARL6IP1", "CENPE", "CITED2", "H2AFZ", "KNSTRN", "STMN1", "ZWINT", "ECT2", "CCNB2", "DYNLL1", "BIRC5", "PTTG1", "UBE2S"

signatures <- list(
  T47_tam_down_signature = c("MT-ATP6", "FMN1", "NDFIP1", "SRP9", "IL6ST", "MRFAP1", "NDUFC2", "TFRC", "FKBP4", "H3F3A", "THBS1", "RPS10", "PPM1K", "ATP5MC3", "MARCKSL1", "MYBL1", "SFPQ", "OLFM1", "HNRNPH3", "PSENEN", "PSMA2", "PCLO", "MINOS1", "B4GALT1", "KRT19", "TIMMDC1", "ARPC5", "PMVK", "ARL3", "MRPL13", "SDF2", "NDUFA8", "CDV3", "MRPL22", "RUNX1", "IGFBP4", "TSFM", "SDSL", "TPBG", "MT-CO3", "GABARAPL2", "GNG11", "LINC00992", "KRT15", "ERH", "TTC14", "AHSA1", "SERPINA6", "REEP5", "LDHA", "SGK3", "MT-CYB", "SLTM", "CA12", "DNAJA1", "TRAPPC4", "GREB1", "HNRNPM", "NDUFB8", "FKBP5", "EIF4A3", "CCT3", "SRSF7", "PFKFB3", "GOLM1", "MAGED2", "PRDX4", "TRPS1", "TMBIM4", "MT-ND2", "HNRNPA2B1", "PGR", "CA8", "XBP1", "SNRPB", "PRLR", "RBBP7", "VDAC3", "CCT5", "PDCD4", "HSP90AA1", "NMU", "DEGS2", "BNIP3", "UBB", "SLC40A1", "JPT1", "TUBA1A", "SPA17", "CALM3", "ENO1", "TUBB4B", "GATA3", "CALM2", "STC2", "H2AFV", "IGFBP5", "ARL6IP1", "CENPE", "CITED2", "H2AFZ", "KNSTRN", "STMN1", "ZWINT", "ECT2", "CCNB2", "DYNLL1", "BIRC5", "PTTG1", "UBE2S")
)

Tumor_cells_53FBAL_NT_Vand <- AddModuleScore_UCell(Tumor_cells_53FBAL_NT_Vand, 
                                                   features=signatures, name=NULL)


head(Tumor_cells_53FBAL_NT_Vand)


# Extract the metadata including the signature scores and orig.ident
signature_data <- FetchData(Tumor_cells_53FBAL_NT_Vand, vars = c("T47_tam_down_signature",
                                                                 "max_factor"))

head(signature_data)


# Subset data for factor 1 and factor 6
factor1_data <- signature_data[signature_data$max_factor == 1, "T47_tam_down_signature"]
factor6_data <- signature_data[signature_data$max_factor == 6, "T47_tam_down_signature"]

# Calculate mean and median for factor 1
mean_factor1 <- mean(factor1_data)
median_factor1 <- median(factor1_data)

# Calculate mean and median for factor 6
mean_factor6 <- mean(factor6_data)
median_factor6 <- median(factor6_data)

# Perform Wilcoxon test between factor 1 and factor 6
wilcox_test <- wilcox.test(factor1_data, factor6_data)

# Display results
mean_factor1
median_factor1
mean_factor6
median_factor6
wilcox_test

# Create violin plots for the three signature scores stored in the metadata
VlnPlot(Tumor_cells_53FBAL_NT_Vand, 
        features = "T47_tam_down_signature", 
        group.by = "max_factor", 
        pt.size = 0)


########
library(Seurat)

DimPlot(Tumor_cells_53FBAL_NT_Vand, group.by = "orig.ident")

Idents(Tumor_cells_53FBAL_NT_Vand) <- "max_factor"
  
# Perform DE analysis using FindMarkers
de_results_NMF1_vs_NMF6 <- FindMarkers(Tumor_cells_53FBAL_NT_Vand, ident.1 = "1", ident.2 = "6", group.by = "max_factor")

write.csv(de_results_NMF1_vs_NMF6, "de_results_NMF1_vs_NMF6.csv")


head()
