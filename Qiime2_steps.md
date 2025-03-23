![image](https://github.com/user-attachments/assets/a7a09bb5-4aa6-40a4-bf86-89026c40b77a)![image](https://github.com/user-attachments/assets/1f70b81e-939d-4b48-a45f-76dcef3fc4d4)# 1. Importing raw data into Qiime2
## Generate manifest.csv
Sequence data are paired end in the format of FASTA with good quality score; therefore, in qiime2 the type will be "SampleData[PairedEndSequencesWithQuality]" and their imput format asigned as PairedEndFastqManifestPhred33.
Before importing, [manifest.csv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/manifest.csv) file must be prepared.

 ```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path short_reads_demux.qza \
  --input-format PairedEndFastqManifestPhred33
```
This might take a while before you get the results. The output of this, is a .qza file that you have already specified it in the command, in this example demuxed-dss.qza. You can then create a visualized file from this artifact, with the following command:

```bash
qiime demux summarize \
  --i-data short_reads_demux.qza \
  --o-visualization short_reads_demux.qzv
```
This [short_reads_demux.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/short_reads_demux.qzv) is a visualized format of short_reads_demux.qza. which you can view it on [qiime2 viewer](https://view.qiime2.org/)). Once you are there you can either drag-and-drop the artifact into the designated area or simpley copy the link to the artifact from this repository and paste it in the box file from the web. Once there, you must come across the following picture:
![image](https://github.com/user-attachments/assets/dbfb9bb8-4dfc-4676-a680-1774ead4cbce)
**Figure 1. Demultiplexed pairedEnd read**    
On this overview page you can see counts of demultiplexed sequences for the entire samples for both forward and reverse reads, with min, median, mean and max and total counts.
![image](https://github.com/user-attachments/assets/352add1e-abed-404b-83dc-12a7d0200684)
**Figure 2. Interacvive plot for demultiplexed pairedEnd reads**    
Understanding this plot is crucial for the denoising step, as it allows you to determine the appropriate truncation length for both forward and reverse reads, ensuring that at least 50% of the reads maintain a quality score (Q) of ≥ 30. You can observe these changes by hovering over the interactive box plots. In this case, the quality of both forward and reverse reads remained consistently high, indicating that minimal (at 246nt) or no truncation may be necessary.

# 2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2.
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs short_reads_demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table dada2/table.qza \
  --o-representative-sequences dada2/rep-seqs.qza \
  --o-denoising-stats dada2/denoising-stats.qza
```
You can convert the denoising-stats.qza file into a [denoising-stats.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/denoising-stats.qzv) file and visualize it on qiime viewer as explained earlier
```bash
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```
![image](https://github.com/user-attachments/assets/11c8e699-d193-46e3-840f-bd042194580f)

**Figure 3. The denoising status of the reads for each sample.**    
You can see the number of filtered reads and also the percentage of non-chimeric sequences after denoising.
The filtered reads and also the percentage of non-chimeric sequences are quite low. You may need to adjust the chimera filtering process in DADA2 or apply an alternative approach as below:
**Perform de novo chimera filtering using VSEARCH in QIIME 2**
```bash
qiime vsearch uchime-denovo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --o-chimeras chimeras.qza \
  --o-nonchimeras rep-seqs-no-chimera.qza \
  --o-stats chimera-stats.qza
```
**Check the chimera filtering summary**
```bash
qiime metadata tabulate \
  --m-input-file chimera-stats.qza \
  --o-visualization chimera-stats.qzv
```
Then use the non-chimeric sequences for further analysis.    
**Filter the feature table to remove chimeric sequences**
```bash
qiime feature-table filter-features \
  --i-table table.qza \
  --m-metadata-file rep-seqs-no-chimera.qza \
  --o-filtered-table filtered-table.qza
```
**Visualize feature table after remove chimeric sequences**
```bash
qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization table-no-chimera.qzv
```
If you drag and drop the [table-no-chimera.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/table-no-chimera.qzv) file in qiime2 view, you can see three main menues; Overview, Interactive Sample Detail and Feature Detail. If you click on Feature detail Detail you can see a slider to the left of the picture which could be changed, based which you can arbiterarily decide, to which depth of reading you can do your rarefaction.
![image](https://github.com/user-attachments/assets/5e1cdc75-645f-4d19-8e0d-e63f5fe5e188)

**Figure 5. ASV table indicating number of samples per treatment and number of ASVs per sample**

# 3. Training a full-length 16S rRNA classifier for taxonomic classification using the Naïve Bayes method in QIIME 2
For taxonomic classifications, you need to have a classifier to which you blast your sequences against to find out which taxonomic groups each sequence belongs to. This is also called reference phylogeny, which is a cruitial step in identifying the marker genes (in this case 16S rRNA) taken from different environmental a in saco samples. In order to do so, there are different 16S rRNA databases, of which Greengens and SILVA are well-known databases for the full length of 16S rRNA genes. You can always download the pre-trained classifiers at the Data Resources of qiime2 website with the follwoing command:
```bash
wget https://data.qiime2.org/2023.7/common/silva-138-99-seqs.qza -O silva_data/silva-138-99-seqs.qza
wget https://data.qiime2.org/2023.7/common/silva-138-99-tax.qza -O silva_data/silva-138-99-tax.qza
```
**Train the downloaded database**
```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva_data/silva-138-99-seqs.qza \
  --i-reference-taxonomy silva_data/silva-138-99-tax.qza \
  --o-classifier silva_data/silva-classifier.qza
```
After you got **silva-classifier.qza** classifier file, you can use it for your taxonomic classifications as follows:
## Taxonomic classification with a confidence score threshold of 90%
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva_data/silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs-no-chimera.qza \
  --p-confidence 0.9 \
  --o-classification taxonomy.qza
```
**Generate taxa bar plot**
```bash
qiime taxa barplot \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-barplot.qzv
```
![image](https://github.com/user-attachments/assets/87310d82-f336-44cc-9638-9a41d90ed04c)
**Figure 6. Taxonomy classification bar plot at genus level**

You can see in the [taxa barplot](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/taxa-barplot.qzv) that most samples have similar microbial compositions, except for one **control sample** and one **RP-20 sample**, which show abnormal patterns. This could be due to low sequencing depth, leading to an inaccurate representation of microbial diversity in these samples. Further analysis is needed to determine whether the microbial compositions in these samples, particularly in the control sample and RP-20 sample, differ significantly in a statistically meaningful way.

To investigate further, we will examine the summary statistics after chimera filtering and perform rarefaction curve analysis in the next steps.    

# 4. Creating a phylogenetic tree using align-to-tree-MAFFT-FastTree
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-no-chimera.qza \
  --o-alignment aligned-rep-seqs-no-chimera.qza \
  --o-masked-alignment masked-aligned-rep-seqs-no-chimera.qza \
  --o-tree unrooted-tree-no-chimera.qza \
  --o-rooted-tree tree-no-chimera.qza
```
**Visualization**
```bash
qiime tools export \
  --input-path tree-no-chimera.qza \
  --output-path exported_tree
```
The Newick file will be in the exported_tree/tree.nwk folder. You can upload [tree.nwk](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/tree.nwk) to [iTOL](https://itol.embl.de/upload.cgi) to view it.
![image](https://github.com/user-attachments/assets/0c598a16-3795-41cf-a63b-3b19e65028e4)


# 5. Rarefraction curve analysis
Rarefaction curves assess sequencing depth sufficiency and microbial diversity saturation in samples. This helps (1) Ensure adequate sequencing depth for reliable diversity estimation (2) Compare samples to detect under-sequenced ones; (3) Evaluate data stability—a plateauing curve indicates sufficient sampling.
**Generate rarefraction curve chart**
```bash
qiime diversity alpha-rarefaction \
  --i-table filtered-table.qza \
  --i-phylogeny tree-no-chimera.qza \
  --p-max-depth 25000 \
  --p-steps 20 \
  --p-iterations 10 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```
You can view the chart [alpha-rarefaction.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/alpha-rarefaction.qzv) by qiime view as method explained previously.
![image](https://github.com/user-attachments/assets/17745783-82e3-4938-aecb-8c9540d007e2)

At a depth of 20000, the rarefaction curve reaches saturation, indicating that increasing reads will not detect many new ASVs, and 15 samples are removed, it means that 100% of the samples are retained. This suggests that normalization at this depth would keep all of the dataset for analysis.    
**Proceed with normalization at a subsampling depth of 20000**
```bash
qiime diversity core-metrics-phylogenetic \
  --i-table filtered-table.qza \
  --i-phylogeny tree-no-chimera.qza \
  --p-sampling-depth 20000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results-20000
```
# 6. Alpha diversity analysis
## Chao1
Chao1 alpha diversity analysis estimates species richness in a sample. It focuses on rare ASVs (appearing only once or twice), helping to assess how many species might be undetected due to limited sequencing depth. This is useful for comparing microbial richness across sample groups.
```bash
qiime diversity alpha \
  --i-table core-metrics-results-20000/rarefied_table.qza \
  --p-metric chao1 \
  --o-alpha-diversity core-metrics-results-20000/chao1_vector.qza
```
```bash
qiime metadata tabulate \
  --m-input-file core-metrics-results-20000/chao1_vector.qza \
  --o-visualization core-metrics-results-20000/chao1_vector.qzv
```
[chao1_vector.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/chao1_vector.qzv) is generated.
![image](https://github.com/user-attachments/assets/3efb70bb-01fb-4f2f-92ec-d1cd1c9d2b30)
Sample 4 has a Chao1 index of 209.25, and Sample 5 has a Chao1 index of 210.86, indicating a lower ASV richness compared to other samples. A statistical test is needed to determine whether this difference is statistically significant among groups.

Analysis Chao1 between groups of treatment
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-20000/chao1_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results-20000/chao1_group_significance.qzv
```
[chao1-group-significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/chao1_group_significance.qzv) is generated.
![image](https://github.com/user-attachments/assets/60b5ead9-2414-433d-acd8-bb15ba935eca)

The Kruskal-Wallis test was performed to compare Chao1 richness across different treatment groups. The overall test statistic H=3.23H=3.23 with a p-value of 0.5196, indicating that there is no statistically significant difference in Chao1 richness among the groups.

Pairwise comparisons using the Kruskal-Wallis test show that all p-values are greater than 0.05, suggesting no significant differences between any two groups. The lowest p-values (e.g., 0.1266) are still far from the significance threshold, further supporting that the observed variations in OTU richness are likely due to random fluctuations rather than meaningful biological differences.

**Conclusion: The results indicate that there is no statistically significant difference in Chao1 richness among the Treatment groups.**

## Shannon index
The Shannon index measures alpha diversity, accounting for both species richness (number of species) and evenness (distribution of species abundances).
    Higher Shannon index → More diverse and evenly distributed microbial community.
    Lower Shannon index → A community dominated by a few species, indicating lower diversity.
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-20000/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results-20000/shannon_group_significance.qzv
```

[shannon-group-significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/shannon_group_significance.qzv) file is generated.
![image](https://github.com/user-attachments/assets/a7a975f0-0a76-4c85-bc84-be0db50dcb18)

Overall: p-value =  	0.597 → There is no significant difference in the Shannon index between groups overall. This means that the microbial diversity structure across all groups is similar.
Pairwise: All p-values > 0.05 → There is no significant difference in the Shannon index between any pair of groups. This indicates that no group has significantly higher or lower microbial diversity compared to the others.    
**Conclusion: The groups have equivalent microbial diversity, suggesting that the grouping factor (e.g., experimental condition) does not strongly influence gut microbiome diversity in this dataset**    

## Pielou's Evenness Index
Pielou's Evenness Index measures the evenness of species distribution in a community. It indicates how evenly the species are distributed, with values ranging from 0 (completely uneven) to 1 (completely even).
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-20000/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results-20000/evenness_group_significance.qzv
```
[evenness_group_significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/evenness_group_significance.qzv) file is generated
![image](https://github.com/user-attachments/assets/28de76cb-ba44-402a-a521-b154e002d13a)

Overall: H =	2.933, p-value = 0.569 → No significant difference in evenness between groups.
Pairwise: All pairwise comparisons have p-value > 0.05, indicating no significant differences in evenness between any pair of groups.    
**Conclusion:There is no significant difference in Pielou's Evenness Index across groups, suggesting that the evenness of microbial distribution is similar in all groups**

## Faith's Phylogenetic Diversity (Faith's PD) 
Faith's PD measures the total branch length of a phylogenetic tree that connects all species in a sample. It reflects both species richness and phylogenetic diversity, considering evolutionary relationships.
A higher Faith’s PD indicates a more diverse microbial community with greater evolutionary variety, while a lower Faith’s PD suggests a more phylogenetically constrained community.
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-20000/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results-20000/faith_pd_group_significance.qzv
```
[faith-pd-group-significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/faith_pd_group_significance.qzv) is generated
![image](https://github.com/user-attachments/assets/aa4215e4-6fa9-4c80-8892-d333ea0ed3d3)

Overall (All groups): H-statistic = 5.566, p-value = 0.2339 → No significant difference in Faith's PD between the groups.    
Pairwise (Between groups):All pairwise comparisons have p-value > 0.05, indicating no significant differences in Faith's PD between any of the groups.    
**Conclusion: There is no significant difference in Faith's PD between the groups, suggesting that the phylogenetic diversity is similar across all groups in this study.**

# 7. Beta diversity analysis
## 1. Bray-Curtis Index
This index is used to measure the difference between two microbial communities. The value ranges from 0 to 1:
    0 means the two communities are completely identical.
    1 means the two communities are completely different.
This index helps us understand the level of difference in microbial species between samples.
**PCoA Plot**
```bash
qiime diversity pcoa \
  --i-distance-matrix core-metrics-results-20000/bray_curtis_distance_matrix.qza \
  --o-pcoa core-metrics-results-20000/bray_curtis_pcoa_results.qza
```
```bash
qiime emperor plot \
  --i-pcoa core-metrics-results-20000/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results-20000/bray_curtis_emperor.qzv
```
![image](https://github.com/user-attachments/assets/7c95322d-b373-45b6-9841-dca732f2cadb)


The PCoA plot [bray_curtis_emperor.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/bray_curtis_emperor.qzv) does not show clear clustering between the treatments,further PERMANOVA analysis for a more detailed examination and reveals significant differences in microbial communities across the groups.
**PERMANOVA analysis**
```bash
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-20000/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --o-visualization core-metrics-results-20000/bray_curtis_significance.qzv \
  --p-method permanova \
  --p-pairwise
```
[bray_curtis_group_significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/bray_curtis_significance.qzv) file is generated.
![image](https://github.com/user-attachments/assets/ac60bc70-1fbd-469f-9256-8dcaad943654)
![image](https://github.com/user-attachments/assets/b74d8b0d-661c-4f52-8eda-c5ddea82cc02)

The overall PERMANOVA test (Bray-Curtis dissimilarity) shows a significant difference between groups (p = 0.017). However, all pairwise comparisons have p-values greater than 0.05, suggesting that while there is a global difference among groups, no specific pairwise comparison shows a statistically significant difference. This could be due to limited sample size or variability within groups.

## Jaccard index
Jaccard index is a measure of similarity between two sets. In microbiome studies, it is used to compare the presence or absence of species across different samples.
    Range: 0 to 1.
        0 means the samples have no species in common (completely different).
        1 means the samples have identical species (completely the same).
Jaccard index focuses on the presence/absence of species, rather than their abundance, making it useful for studies where the presence of specific species is more important than their relative abundance.
```bash
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-20000/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --o-visualization core-metrics-results-20000/jaccard_pairwise_significance.qzv \
  --p-method permanova \
  --p-pairwise
```
[jaccard_pairwise_significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/jaccard_pairwise_significance.qzv) is generated
![image](https://github.com/user-attachments/assets/4bd9a800-afb9-49c5-860f-2a578d62bc24)
![image](https://github.com/user-attachments/assets/0c172bff-caa1-4707-93c7-fdb4fdfa44c0)
    p-value = 0.13: Since the p-value is greater than 0.05, this indicates that there is no significant difference in the microbial community structure between the Treatment groups.
    Test Statistic (pseudo-F) = 1.006371: This is the statistical value used to assess the difference between groups. However, the result is not strong enough to indicate a clear distinction, as reflected by the large p-value.    
**Conclusion:There is no significant differentiation between the Treatment groups in this dataset based on the Jaccard index and PERMANOVA analysis.**

## Weighted UniFrac 
Weighted UniFrac is better at detecting differences in community structure by considering the relative abundance of each species, not just whether or not they are present in a sample.
```bash
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-20000/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --o-visualization core-metrics-results-20000/weighted_unifrac_significance.qzv \
  --p-method permanova \
  --p-pairwise
```
[weighted_unifrac_significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/weighted_unifrac_significance.qzv) file is generated
![image](https://github.com/user-attachments/assets/9e807173-6742-4ae7-a5d9-2af1bd28e229)
![image](https://github.com/user-attachments/assets/445616c1-eb2b-45fe-801d-f59fa51f64ea)

    p-value = 0.113: The p-value is greater than 0.05, indicating that there is no significant difference in the microbial community structure between the treatment groups.
    Test Statistic (pseudo-F) = 1.842749: suggests some variation among groups, but it is not strong enough to reach statistical significance. This could be due to small sample size or high within-group variability.
**Conclusion: There is no significant difference in the microbial community structure between the treatment groups based on the Weighted UniFrac distance matrix and the PERMANOVA test.**

## Unweighted UniFrac
Unweighted UniFrac measures the differences between microbial communities based on the presence or absence of species, without considering their abundance or frequency. It is useful for comparing the community structure and detecting shifts in species diversity between samples, especially when the focus is not on the abundance of each species.
```bash
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-20000/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --o-visualization core-metrics-results-20000/unweighted_unifrac_significance.qzv \
  --p-method permanova \
  --p-pairwise
```
[unweighted_unifrac_significance.qzv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Qiime_steps/unweighted_unifrac_significance.qzv) file is generated.
![image](https://github.com/user-attachments/assets/1d9d6602-5d73-4701-9836-b8f34525addf)
![image](https://github.com/user-attachments/assets/3ae0a9ca-d229-4ba5-98c2-3b594eede721)


The results from unweighted_unifrac_significance.qzv (test statistic = 0.91104, p-value = 0.901) indicate that there is no significant difference in microbial community composition between the treatment groups based on unweighted UniFrac distance. The high p-value suggests that any observed variations are likely due to random chance rather than a meaningful biological difference.
**Conclusion: There is no significant difference in microbial community structure across the treatment groups based on the unweighted UniFrac distance, suggesting that the presence/absence of species is not strongly influenced by the treatments in this dataset.**

# 8. Export file from Qimme2 for R steps
```bash
qiime tools export \
  --input-path table.qza \
  --output-path exported_feature_table
```
The above command will export the table.qza file as a temporary file (like a .biom file) in the exported_feature_table directory.

Then, to convert the .biom file into a data table format that R can handle (CSV/TSV) - [feature_table.tsv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/feature-table.tsv), you can use the following command:
```bash
biom convert \
  -i exported_feature_table/feature-table.biom \
  -o feature_table.csv \
  --to-tsv
```

```bash
qiime taxa export \
  --i-classification taxonomy.qza \
  --output-dir taxonomy_exported
```
The above command will export the [taxanomy.tsv](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/R_steps/taxonomy.tsv) file in the taxonomy_exported directory.
