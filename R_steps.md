# 1. Differential abundance analysis of taxa by ANCOMBC
## 1.1 Set up working environment
```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
```

## 1.2. Dataset
```bash
feature_table <- read.table("feature-table.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
```

## 1.3. Converting Data to Phyloseq Object
```bash
library(phyloseq)

# Create an ASV table object
asv_phyloseq <- asv_table(feature_table, taxa_are_rows = TRUE)

# Create a sample_data object from metadata
sample_metadata <- sample_data(metadata)

# Check the structure of the ASV table
str(asv_phyloseq)

# Combine ASV table and metadata into a Phyloseq object
ps <- phyloseq(asv_phyloseq, sample_metadata)
```

## 1.4. Matching Sample Names
```bash
# Check sample names in ASV table and metadata
sample_names(asv_phyloseq)
rownames(sample_metadata)

# Rename sample names in ASV table to match metadata
sample_names(asv_phyloseq) <- rownames(sample_metadata)

# Verify the sample names in the ASV table
sample_names(asv_phyloseq)
```
[pairwise_comparisons_results.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/pairwise_comparisons_results.csv) file is created.

## 1.5.  ANCOM-BC2 Analysis
```bash
res_ancombc2 <- ancombc2(
  data = ps,
  fix_formula = "Treatment",  # The variable for analysis
  p_adj_method = "holm",      # Method for adjusting p-values
  alpha = 0.05,               # Significance level
  pairwise = TRUE,            # Enable pairwise comparisons
  group = "Treatment"         # The group variable
)
```

## 1.6 Add a pseudo-count to the ASV table to handle zero values
```bash
# Add a pseudo-count to the ASV table to handle zero values
asv_table_with_pseudo <- asv_table(ps) + 1

# Create a new phyloseq object with the updated OTU table
ps_with_pseudo <- phyloseq(asv_table(asv_table_with_pseudo, taxa_are_rows = TRUE), sample_data(ps))
```

## 1.7 Save the ANCOM-BC2 results to a CSV file
```bash
ancombc2_results <- res_ancombc2$res
write.csv(ancombc2_results, "ancombc2_results.csv", row.names = TRUE)

# Save pairwise comparison results to a CSV file
write.csv(res_ancombc2$res_pair, "pairwise_comparisons_results.csv", row.names = TRUE)
```
[pairwise_comparisons_results.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/pairwise_comparisons_results.csv) file is created.

# 2. Analysis of biomarkers for gut microbiome
## 2.1. Filtering Taxa with Significant p-value<0.05
```bash
# Filter taxa with p-value < 0.05
df_filtered <- df_long_lfc_pvalue %>%
    filter(pvalue < 0.05)

# Save filtered results to a CSV file
write.csv(df_filtered, "filtered_taxon_results.csv", row.names = TRUE)
```
## 2.2. Visualization of Log Fold Changes of Taxa with Significant p-value<0.05
```bash
library(dplyr)
library(ggplot2)

# Assuming merged_data has been processed and contains columns: treatment, lfc, final_taxon
# Filter out rows where final_taxon is "uncultured", "unassigned", "uncultured_soil", "uncultured_rumen"
merged_data_filtered <- merged_data %>% 
  filter(!(tolower(final_taxon) %in% c("uncultured", "unassigned", "uncultured_soil", "uncultured_rumen")))

# Convert treatment column to a factor with the desired order
merged_data_filtered$treatment <- factor(merged_data_filtered$treatment,
                                         levels = c("RP.5", "RP.10", "RP.20", "RP.40", "RP.X"))

# If lfc is initially in natural log fold change, convert it to log2 fold change
merged_data_filtered <- merged_data_filtered %>% mutate(log2fc = lfc / log(2))

# Create Volcano Plot:
# - x-axis: log2 fold change
# - y-axis: final_taxon
# - Color and shape points by treatment
p <- ggplot(merged_data_filtered, aes(x = log2fc, y = final_taxon, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("RP.5" = "red",
                                "RP.10" = "blue",
                                "RP.20" = "green",
                                "RP.40" = "orange",
                                "RP.X" = "purple")) +
  scale_shape_manual(values = c(16, 17, 15, 18, 3)) +
  labs(x = "Log2 Fold Change (pvalue < 0.05)",
       y = "Final Taxon",
       color = "Treatment",
       shape = "Treatment") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

print(p)
```
```bash
# Save the plot as a PNG file
ggsave("lfc_by_treatment.png", plot = last_plot(), width = 8, height = 6)
```
[lfc_by_treatment.png](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/lfc_by_treatment.png.png) is a visualization of taxa that have significant log fold change values.
![image](https://github.com/user-attachments/assets/642860ba-d86d-4151-b122-a32078833e0b)

