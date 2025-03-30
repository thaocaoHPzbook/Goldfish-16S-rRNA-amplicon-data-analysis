# 1. Differential abundance analysis of taxa by ANCOMBC
## 1.1 Set up working environment
```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
```

## 1.2. Read and Prepare Data
```bash
feature_table <- read.table("feature-table.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
```

## 1.3. Convert Data into Phyloseq Object
```bash
library(phyloseq)

# Create an ASV table object
otu_phyloseq <- otu_table(feature_table, taxa_are_rows = TRUE)

# Create a sample_data object from metadata
sample_metadata <- sample_data(metadata)

# Check the structure of the ASV table
str(otu_phyloseq)

# Combine ASV table and metadata into a Phyloseq object
ps <- phyloseq(otu_phyloseq, sample_metadata)
```

## 1.4. Matching Sample Names
```bash
# Check sample names in the OTU table
sample_names(otu_phyloseq)

# Check sample names in the metadata
rownames(sample_metadata)

# Rename OTU table sample names to match metadata sample names
sample_names(otu_phyloseq) <- rownames(sample_metadata)

# Removes any X prefixes that might have been introduced when importing data
sample_names(otu_phyloseq) <- gsub("^X", "", sample_names(otu_phyloseq))

# Ensure sample names are numeric (if applicable)
sample_names(otu_phyloseq) <- as.character(as.numeric(sample_names(otu_phyloseq)))

# Reorder samples in the OTU table to match metadata order
otu_phyloseq <- prune_samples(rownames(sample_metadata), otu_phyloseq)

# Transform sample counts without changing OTU values
otu_phyloseq <- transform_sample_counts(otu_phyloseq, identity)

# Reapply correct sample names
sample_names(otu_phyloseq) <- rownames(sample_metadata)

# Check sample names after processing
print(sample_names(otu_phyloseq))
print(rownames(sample_metadata))

# Create a phyloseq object combining OTU table and metadata
ps <- phyloseq(otu_phyloseq, sample_metadata)
```

## 1.5 Add a pseudo-count to the OUT table to handle zero values
```bash
# Add a pseudo-count to the ASV table to handle zero values
otu_table_with_pseudo <- otu_table(ps) + 0.5

# Create a new phyloseq object with pseudo-count added OTU table
ps_with_pseudo <- phyloseq(otu_table(otu_table_with_pseudo, taxa_are_rows = TRUE), sample_data(ps))

# Check the updated phyloseq object
ps_with_pseudo
```
## 1.6 Convert Phyloseq Object into OTU Table and Metadata
**Extract OTU table and metadata from phyloseq object**
```bash
otu_data <- otu_table(ps_with_pseudo)
sample_data <- sample_data(ps_with_pseudo)
```

##  1.7 Run ANCOM-BC2 with Pseudo-Count Adjusted Data
```bash
res_ancombc2 <- ancombc2(
  data = ps_with_pseudo,    # Phyloseq object with pseudo-counts
  fix_formula = "Treatment", # Variable for comparison
  p_adj_method = "holm",     # Multiple testing correction
  alpha = 0.05,              # Significance level
  pairwise = TRUE,           # Enable pairwise comparison
  group = "Treatment"        # Grouping variable
)
```

## 1.8 Save the ANCOM-BC2 results to a CSV file
```bash
ancombc2_results <- res_ancombc2$res
write.csv(ancombc2_results, "ancombc2_results.csv", row.names = TRUE)

# Save pairwise comparison results to a CSV file
write.csv(res_ancombc2$res_pair, "pairwise_comparisons_results.csv", row.names = TRUE)
```
[ancombc2_results.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/ancombc2_results.csv), [pairwise_comparisons_results.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/pairwise_comparisons_results.csv) files are created.

# 2. Analysis of biomarkers for gut microbiome
## 2.1. Filtering Taxa with Significant p-value<0.05
```bash
# Filter taxa with p-value < 0.05, q-value<0.5, diff=TRUE, and passed_ss=TRUE
library(tidyr)
library(dplyr)
filtered_ancombc2_results <- ancombc2_results %>%
  filter(
    (p_TreatmentRP.5 < 0.05 & q_TreatmentRP.5 < 0.05 & diff_TreatmentRP.5 == TRUE & passed_ss_TreatmentRP.5 == TRUE) |
    (p_TreatmentRP.10 < 0.05 & q_TreatmentRP.10 < 0.05 & diff_TreatmentRP.10 == TRUE & passed_ss_TreatmentRP.10 == TRUE) |
    (p_TreatmentRP.20 < 0.05 & q_TreatmentRP.20 < 0.05 & diff_TreatmentRP.20 == TRUE & passed_ss_TreatmentRP.20 == TRUE) |
    (p_TreatmentRP.40 < 0.05 & q_TreatmentRP.40 < 0.05 & diff_TreatmentRP.40 == TRUE & passed_ss_TreatmentRP.40 == TRUE)
  )
# Save filtered results to a CSV file
write.csv(df_filtered, "filtered_ancombc2_results.csv", row.names = TRUE)
```

**Merge taxomomy.tsv with filtered_ancombc2_results.csv**
```bash
library(dplyr)
library(readr)

# Read taxonomy.tsv file
taxonomy <- read_tsv("taxonomy.tsv")

# Read filtered_ancombc2_parirwise_results.csv file
filtered_results <- read_csv("filtered_ancombc2_results.csv")

# Join data
merged_results <- filtered_results %>%
  left_join(taxonomy, by = c("taxon" = "Feature ID"))

# Check the results
head(merged_results)

# Replace taxon with Taxon name
final_results <- merged_results %>%
  select(-taxon) %>%
  rename(taxon = Taxon)

# Exporting csv file
write_csv(final_results, "filtered_ancombc2_results_with_taxonomy.csv")
```
[filtered_ancombc2_results_with_taxonomy.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/filtered_ancombc2_results_with_taxonomy.csv) is generated.

## 2.2. Visualization of Log Fold Changes of Taxa with Significant p-value<0.05
```bash
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(readr)

# Đọc dữ liệu
filtered_ancombc2_results <- read_csv("filtered_ancombc2_results_with_taxonomy.csv")

# Trích xuất genus từ taxon
filtered_ancombc2_results$genus <- str_extract(filtered_ancombc2_results$taxon, "g__[^;]+")

# Chuyển dữ liệu từ wide format sang long format
filtered_long <- filtered_ancombc2_results %>%
  pivot_longer(cols = starts_with("lfc_TreatmentRP"),
               names_to = "treatment",
               values_to = "lfc") %>%
  mutate(treatment = str_replace(treatment, "lfc_Treatment", "")) %>%
  pivot_longer(cols = starts_with("p_TreatmentRP"),
               names_to = "p_treatment",
               values_to = "p_value") %>%
  filter(str_replace(p_treatment, "p_Treatment", "") == treatment) %>% # Ghép đúng lfc với p-value
  select(-p_treatment)

# Xác định dấu * cho các điểm có p-value < 0.05
filtered_long$significance <- ifelse(filtered_long$p_value < 0.05, "*", "")

# Vẽ biểu đồ
p <- ggplot(filtered_long, aes(x = genus, y = lfc, fill = treatment)) +
    geom_col(position = position_dodge(width = 0.7),  # Xếp cột sát nhau
             width = 0.6,  # Điều chỉnh kích thước cột
             color = "black") +  
    geom_text(aes(label = significance, y = lfc + sign(lfc) * 0.4),  # Đẩy dấu * cao hơn
              position = position_dodge(width = 0.7), 
              size = 6, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +  # Đường tham chiếu
    scale_fill_manual(values = c("RP.5" = "red",
                                 "RP.10" = "blue",
                                 "RP.20" = "green",
                                 "RP.40" = "orange",
                                 "RP.X" = "purple")) +
    labs(x = "Genus",
         y = "Log2 Fold Change",
         fill = "Treatment") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Xoay label trục X
          legend.position = "right")  # Giữ legend bên phải

# Hiển thị biểu đồ
print(p)
```
```bash
# Save the plot as a PNG file
ggsave("lfc_by_treatment.png", plot = last_plot(), width = 8, height = 6)
```
[lfc_by_treatment.png](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/lfc_by_treatment.png) is a visualization of taxa that have significant log fold change values.
![image](https://github.com/user-attachments/assets/ff26fbfd-5206-4a6e-ab86-1954ac640df5)



## 2.3  Visualization of Log Fold Changes of Taxa with Significant p-value<0.05 - pairwise comparision
```bash
# Filter taxa with p-value < 0.05, q-value<0.5, diff=TRUE, and passed_ss=TRUE
filtered_ancombc2_results <- ancombc2_results %>%
  filter(
    (.data[["p_TreatmentRP-20_TreatmentRP-10"]] < 0.05 &
     .data[["q_TreatmentRP-20_TreatmentRP-10"]] < 0.05 &
     .data[["diff_TreatmentRP-20_TreatmentRP-10"]] == TRUE &
     .data[["passed_ss_TreatmentRP-20_TreatmentRP-10"]] == TRUE) |
    
    (.data[["p_TreatmentRP-40_TreatmentRP-10"]] < 0.05 &
     .data[["q_TreatmentRP-40_TreatmentRP-10"]] < 0.05 &
     .data[["diff_TreatmentRP-40_TreatmentRP-10"]] == TRUE &
     .data[["passed_ss_TreatmentRP-40_TreatmentRP-10"]] == TRUE) |
    
    (.data[["p_TreatmentRP-5_TreatmentRP-10"]] < 0.05 &
     .data[["q_TreatmentRP-5_TreatmentRP-10"]] < 0.05 &
     .data[["diff_TreatmentRP-5_TreatmentRP-10"]] == TRUE &
     .data[["passed_ss_TreatmentRP-5_TreatmentRP-10"]] == TRUE) |
    
    (.data[["p_TreatmentRP-40_TreatmentRP-20"]] < 0.05 &
     .data[["q_TreatmentRP-40_TreatmentRP-20"]] < 0.05 &
     .data[["diff_TreatmentRP-40_TreatmentRP-20"]] == TRUE &
     .data[["passed_ss_TreatmentRP-40_TreatmentRP-20"]] == TRUE) |
    
    (.data[["p_TreatmentRP-5_TreatmentRP-20"]] < 0.05 &
     .data[["q_TreatmentRP-5_TreatmentRP-20"]] < 0.05 &
     .data[["diff_TreatmentRP-5_TreatmentRP-20"]] == TRUE &
     .data[["passed_ss_TreatmentRP-5_TreatmentRP-20"]] == TRUE) |
    
    (.data[["p_TreatmentRP-5_TreatmentRP-40"]] < 0.05 &
     .data[["q_TreatmentRP-5_TreatmentRP-40"]] < 0.05 &
     .data[["diff_TreatmentRP-5_TreatmentRP-40"]] == TRUE &
     .data[["passed_ss_TreatmentRP-5_TreatmentRP-40"]] == TRUE)
  )

# Save filtered results to a CSV file
write.csv(df_filtered, "filtered_ancombc2_parirwise_results.csv", row.names = TRUE)
```
[filtered_ancombc2_parirwise_results.csv"](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/filtered_ancombc2_pairwise_results_with_taxonomy.csv) is generated.
**Merge taxomomy.tsv with filtered_ancombc2_parirwise_results.csv**
```bash
library(dplyr)
library(readr)

# Read taxonomy.tsv file
taxonomy <- read_tsv("taxonomy.tsv")

# Read filtered_ancombc2_parirwise_results.csv file
filtered_results <- read_csv("filtered_ancombc2_pairwise_results.csv")

# Join data
merged_results <- filtered_results %>%
  left_join(taxonomy, by = c("taxon" = "Feature ID"))

# Check the results
head(merged_results)

# Replace taxon with Taxon name
final_results <- merged_results %>%
  select(-taxon) %>%
  rename(taxon = Taxon)

# Exporting csv file
write_csv(final_results, "filtered_ancombc2_pairwise_results_with_taxonomy.csv")
```


**Plotting**
```bash
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

# 1. Danh sách các so sánh cần thiết
target_cols <- c(
    "lfc_TreatmentRP.20_TreatmentRP.10",
    "lfc_TreatmentRP.40_TreatmentRP.10",
    "lfc_TreatmentRP.5_TreatmentRP.10",
    "lfc_TreatmentRP.40_TreatmentRP.20",
    "lfc_TreatmentRP.5_TreatmentRP.20",
    "lfc_TreatmentRP.5_TreatmentRP.40"
)

pval_cols <- str_replace(target_cols, "lfc_", "p_")  # Tạo danh sách cột p-value

# 2. Lọc taxon đến cấp genus và chuyển sang dạng long format
long_results <- filtered_ancombc2_results %>%
    filter(!is.na(taxon)) %>%
    pivot_longer(
        cols = all_of(target_cols),
        names_to = "comparison",
        values_to = "lfc"
    ) %>%
    mutate(
        comparison = str_replace(comparison, "lfc_", ""),
        genus = str_extract(taxon, "g__[^;]+")
    ) %>%
    filter(!is.na(genus)) %>%  # Loại bỏ taxon chưa phân loại tới genus
    inner_join(  # Dùng inner_join để chỉ lấy các giá trị phù hợp
        filtered_ancombc2_results %>%
            pivot_longer(
                cols = all_of(pval_cols),
                names_to = "comparison",
                values_to = "pvalue"
            ) %>%
            mutate(comparison = str_replace(comparison, "p_", "")),
        by = c("taxon", "comparison")
    ) %>%
    group_by(genus, comparison) %>%  # Nhóm theo genus và comparison
    summarise(
        lfc = mean(lfc, na.rm = TRUE),  # Lấy trung bình nếu có nhiều giá trị lặp
        pvalue = mean(pvalue, na.rm = TRUE)
    ) %>%
    mutate(significant = ifelse(pvalue < 0.05, "*", "")) %>%
    ungroup()

# 3. Vẽ biểu đồ cột
p <- ggplot(long_results, aes(x = genus, y = lfc, fill = comparison)) +
    geom_col(position = position_dodge(width = 0.8), color = "black") +  # Cột không bị trùng
    geom_text(aes(label = significant), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 6, color = "black") +  # Giữ dấu * đúng vị trí
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +  # Đường tham chiếu
    scale_fill_manual(
        values = c(
            "TreatmentRP.20_TreatmentRP.10" = "#E41A1C",
            "TreatmentRP.40_TreatmentRP.10" = "#377EB8",
            "TreatmentRP.5_TreatmentRP.10"  = "#4DAF4A",
            "TreatmentRP.40_TreatmentRP.20" = "#984EA3",
            "TreatmentRP.5_TreatmentRP.20"  = "#FF7F00",
            "TreatmentRP.5_TreatmentRP.40"  = "#A65628"
        )
    ) +
    labs(
        x = "Genus",
        y = "Log2 Fold Change",
        fill = "Comparison"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Xoay label trục X để dễ đọc
        legend.position = "bottom"
    )

# Hiển thị biểu đồ
print(p)
```
**Save the plot**
```bash
# Save the plot as a PNG file
ggsave("lfc_by_treatment_pairwise.png", plot = last_plot(), width = 8, height = 6)
```
[lfc_by_treatment_pairwise.png]
![image](https://github.com/user-attachments/assets/99fb8b7a-45b8-4e55-80fa-a3396530a5a1)




