# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(openxlsx)  # For saving to Excel
library(ggplot2)
library(ggrepel)  # For better label placement

# Load and merge datasets
nutlin_data <- load_mixseq_data("Idasanutlin_24hr_expt1", "Idasanutlin")
dmso_data <- load_mixseq_data("DMSO_24hr_expt1", "DMSO")
sc_data <- merge(nutlin_data, dmso_data)

# Check unique values in condition_TP53_status
unique(sc_data@meta.data$condition_TP53_status)

# Subset sc_data based on condition_TP53_status containing "Idasanutlin"
sc_data_clean <- sc_data[, grepl("Idasanutlin", sc_data@meta.data$condition_TP53_status)]

# Check dimensions or unique values in sc_data_clean
dim(sc_data_clean)  # Check dimensions
unique(sc_data_clean@meta.data$condition_TP53_status)  # Check unique values

# Normalize and process data (if not already done)
npcs <- length(unique(sc_data_clean$singlet_ID)) * 2
sc_data_clean <- NormalizeData(sc_data_clean, verbose = FALSE)
sc_data_clean <- ScaleData(sc_data_clean, verbose = FALSE)
sc_data_clean <- FindVariableFeatures(sc_data_clean, nfeatures = 5000, verbose = FALSE)
sc_data_clean <- RunPCA(sc_data_clean, features = VariableFeatures(sc_data_clean), npcs = npcs, verbose = FALSE)
sc_data_clean <- RunUMAP(sc_data_clean, dims = 1:npcs, n.neighbors = 5, verbose = FALSE)

# Create UMAP plot
unique(sc_data_clean@meta.data$condition_TP53_status)

umap_plot <- DimPlot(sc_data_clean, group.by = "condition_TP53_status", 
                     cols = c("Idasanutlin  TP53_WT" = "red", "Idasanutlin  TP53_MUT" = "black")) +
  labs(title = "UMAP of Cell Lines by Condition and TP53 Status") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")


print(umap_plot)

# Extract UMAP coordinates from the cleaned data
umap_data <- Embeddings(sc_data_clean, "umap")
umap_df <- as.data.frame(umap_data)
umap_df$TP53_status <- sc_data_clean$TP53_status
umap_df$condition <- sc_data_clean$condition

# Save UMAP coordinates to Excel
write.xlsx(umap_df, "umap_coordinates.xlsx", rowNames = FALSE)

# Calculate descriptive statistics for UMAP_1 and UMAP_2
umap_stats <- umap_df %>%
  group_by(TP53_status) %>%
  summarise(
    mean_UMAP1 = mean(umap_1),
    median_UMAP1 = median(umap_1),
    sd_UMAP1 = sd(umap_1),
    mean_UMAP2 = mean(umap_2),
    median_UMAP2 = median(umap_2),
    sd_UMAP2 = sd(umap_2)
  )

print(umap_stats)

# Create a UMAP plot with TP53 status and specified colors
umap_plot <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = TP53_status)) +
  geom_point() +
  scale_color_manual(values = c("TP53_WT" = "red", "TP53_MUT" = "black")) +
  labs(title = "UMAP of Cell Lines by TP53 Status",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")

print(umap_plot)

# Perform statistical tests
# Compute mean and sd for umap_1 and umap_2 by TP53_status
umap_summary <- umap_df %>%
  group_by(TP53_status) %>%
  summarise(
    mean_umap1 = mean(umap_1),
    sd_umap1 = sd(umap_1),
    mean_umap2 = mean(umap_2),
    sd_umap2 = sd(umap_2)
  )

# Print the summary
print(umap_summary)

# Plot mean and sd of UMAP dimensions by TP53 status
umap_summary_plot <- ggplot(umap_summary, aes(x = TP53_status)) +
  geom_errorbar(aes(ymin = mean_umap1 - sd_umap1, ymax = mean_umap1 + sd_umap1), width = 0.2) +
  geom_errorbar(aes(ymin = mean_umap2 - sd_umap2, ymax = mean_umap2 + sd_umap2), width = 0.2) +
  geom_point(aes(y = mean_umap1), color = "red", size = 3) +
  geom_point(aes(y = mean_umap2), color = "black", size = 3) +
  labs(title = "Mean and Standard Deviation of UMAP Dimensions by TP53 Status",
       x = "TP53 Status",
       y = "UMAP Dimensions") +
  theme_minimal()

print(umap_summary_plot)

library(writexl)

# Save to Excel
write_xlsx(umap_summary, "umap_summary.xlsx")
# Save to CSV
write.csv(umap_summary, "umap_summary.csv", row.names = FALSE)


# Create the volcano plot

#Explanation
#Log2 Fold Change Calculation: The log2fc vector calculates the mean expression of each gene in the TP53_WT and TP53_MUT groups and then takes the log2 of their ratio.
#Combining Results: The results_df data frame combines the log2 fold changes and p-values.
#Volcano Plot: The ggplot code creates a volcano plot, highlighting significant genes (adjusted p-value < 0.05) in red.

# Statistical test for expression of top genes between TP53 statuses

# Get the counts data from the correct slot
counts_data <- GetAssayData(sc_data, assay = "RNA", slot = "counts.1")

# Get the cell names from counts_data
counts_cell_names <- colnames(counts_data)

# Identify valid cells by matching cell names
valid_cells <- colnames(sc_data) %in% counts_cell_names

# Debugging steps
print(dim(counts_data))  # Check dimensions of counts_data
print(length(valid_cells))  # Check length of valid_cells
print(class(valid_cells))  # Should be 'logical'
print(length(colnames(sc_data)))  # Number of columns in sc_data
print(sum(valid_cells))  # Number of valid cells
# Example: Assuming you have some way to define top genes (replace with your actual method)
top_genes <- rownames(counts_data)[1:100]  # Selecting top 100 genes for example

# Calculate p-values and log2 fold changes for each gene
gene_stats <- sapply(top_genes, function(gene) {
  wt_cells <- tp53_status_filtered == "TP53_WT"
  mut_cells <- tp53_status_filtered == "TP53_MUT"
  
  print(paste("Gene:", gene))
  print(paste("WT Cells:", sum(wt_cells)))
  print(paste("Mut Cells:", sum(mut_cells)))
  
  if (sum(wt_cells) > 0 && sum(mut_cells) > 0) {
    p_value <- t.test(counts_data[gene, wt_cells], counts_data[gene, mut_cells])$p.value
    log2fc <- log2(mean(counts_data[gene, mut_cells]) / mean(counts_data[gene, wt_cells]))
    c(log2fc, p_value)
  } else {
    c(NA, NA)
  }
})

# Create a data frame with the results, handling NA values
results_df <- data.frame(
  Gene = top_genes,
  log2FoldChange = ifelse(is.na(gene_stats[1, ]), NA, gene_stats[1, ]),
  P_Value = ifelse(is.na(gene_stats[2, ]), NA, gene_stats[2, ]),
  p_adj = ifelse(is.na(gene_stats[2, ]), NA, p.adjust(gene_stats[2, ], method = "BH"))
)

# Filter out rows with NA values
results_df <- results_df[complete.cases(results_df), ]

# Print the filtered results_df
print(results_df)

# Generate the volcano plot with filtered results
library(ggplot2)
library(ggrepel)

volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(P_Value))) + 
  geom_point(aes(color = p_adj < 0.05), alpha = 0.7, size = 3) + 
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("grey", "red")) +
  theme(legend.position = "none") +
  geom_text_repel(data = subset(results_df, p_adj < 0.05), aes(label = Gene), size = 3, max.overlaps = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

print(volcano_plot)


###END##

###Reviewers Comments
# Load necessary libraries (assuming they are already loaded)

# Perform UMAP and statistical test
sc_data_clean <- RunUMAP(sc_data_clean, dims = 1:npcs, n.neighbors = 5, verbose = FALSE)

# Create UMAP plot
umap_plot <- DimPlot(sc_data_clean, group.by = "condition_TP53_status", 
                     cols = c("Idasanutlin  TP53_WT" = "red", "Idasanutlin  TP53_MUT" = "black")) +
  labs(title = "UMAP of Cell Lines by Condition and TP53 Status") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")

print(umap_plot)

# Perform statistical test between TP53_WT and TP53_MUT UMAP coordinates
umap_test <- t.test(umap_df$umap_1 ~ umap_df$TP53_status)

print(umap_test)

# Save UMAP coordinates to Excel
write.xlsx(umap_df, "umap_coordinates.xlsx", rowNames = FALSE)

# Save to CSV
write.csv(umap_summary, "umap_summary.csv", row.names = FALSE)

