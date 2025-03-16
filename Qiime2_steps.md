# 1. Importing raw data into Qiime2
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
This **short_reads_demux.qzv** is a visualized format of short_reads_demux.qza. which you can view it on {qiime2 viewer}(https://view.qiime2.org/)). Once you are there you can either drag-and-drop the artifact into the designated area or simpley copy the link to the artifact from this repository and paste it in the box file from the web. Once there, you must come across the following picture:
![image](https://github.com/user-attachments/assets/dbfb9bb8-4dfc-4676-a680-1774ead4cbce)
**Figure 1. Demultiplexed pairedEnd read**
On this overview page you can see counts of demultiplexed sequences for the entire samples for both forward and reverse reads, with min, median, mean and max and total counts.
![image](https://github.com/user-attachments/assets/352add1e-abed-404b-83dc-12a7d0200684)
**Figure 2. Interacvive plot for demultiplexed pairedEnd reads**
Understanding this plot is crucial for the denoising step, as you need to determine the truncation length for both forward and reverse reads in a way that ensures at least 50% of the reads have a quality score (Q) ≥ 30. You can observe these changes by hovering over the interactive box plots. In this case, the quality of both forward and reverse reads starts to decline significantly after approximately 220 nt.

# 2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2.
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs short_reads_demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```
You can convert the denoising-stats.qza file into a .qzv file and visualize it on qiime viewer as explained earlier
```bash
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```
![image](https://github.com/user-attachments/assets/6de3099b-4481-401a-acf4-d81eeb8ddb72)
Figure 3. The denoising status of the reads for each sample. You can see the number of filtered reads and also the percentage of non-chimeric sequences after denoising.
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
  --o-filtered-table table-no-chimera.qza
```
**Visualize feature table after remove chimeric sequences**
```bash
qiime feature-table summarize \
  --i-table table-no-chimera.qza \
  --o-visualization table-no-chimera.qzv
```
If you drag and drop the **table-no-chimera.qzv** file in qiime2 view, you can see three main menues; Overview, Interactive Sample Detail and Feature Detail. If you click on Feature detail Detail you can see a slider to the left of the picture which could be changed, based which you can arbiterarily decide, to which depth of reading you can do your rarefaction.
![image](https://github.com/user-attachments/assets/5ba5e1a4-054f-4989-b556-cb4f779ee42e)
**Figure 5. ASV table indicating number of samples per treatment and number of ASVs per sample**

# 3. Training a primer-based region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
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
## Taxonomic clasification
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva_data/silva-classifier.qza \
  --i-reads rep-seqs-no-chimera.qza \
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
![image](https://github.com/user-attachments/assets/75829c6a-9148-4f92-88c2-a55da47c00b5)
**Figure 6. Taxonomy classification bar plot at level 6**

You can see in the taxa barplot that most samples have similar microbial compositions, except for one **control sample** and one **RP-20 sample**, which show abnormal patterns. This could be due to low sequencing depth, leading to an inaccurate representation of microbial diversity in these samples.

To investigate further, we will examine the summary statistics after chimera filtering and perform rarefaction curve analysis in the next steps. 

# Creating a phylogenetic tree using align-to-tree-MAFFT-FastTree
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
The Newick file will be in the exported_tree/tree.nwk folder. You can upload it to [iTOL](https://itol.embl.de/upload.cgi) to view it.
![image](https://github.com/user-attachments/assets/20474d10-c2fd-4cf6-a923-b7e20afb4f00)

# Rarefraction curve analysis
Rarefaction curves assess sequencing depth sufficiency and microbial diversity saturation in samples. This helps (1) Ensure adequate sequencing depth for reliable diversity estimation (2) Compare samples to detect under-sequenced ones; (3) Evaluate data stability—a plateauing curve indicates sufficient sampling.
**Generate rarefraction curve chart**
```bash
qiime diversity alpha-rarefaction \
  --i-table filtered-table.qza \
  --i-phylogeny tree-no-chimera.qza \
  --p-max-depth 14007 \
  --p-steps 20 \
  --p-iterations 10 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```
You can view the chart **alpha-rarefaction.qzv** by qiime view as method explained previously.
![image](https://github.com/user-attachments/assets/7bdadcb9-a4c5-4cfd-a2f4-f231d4f37a40)
At a depth of 4000, the rarefaction curve reaches saturation, indicating that increasing reads will not detect many new ASVs.

![image](https://github.com/user-attachments/assets/0efd2bda-dbdf-4a18-b38b-58b9512bfedd)
At a depth of 4000, 2 out of 15 samples are removed, it means that over 80% of the samples (13/15) are retained. This suggests that normalization at this depth would keep most of the dataset for analysis.
**Proceed with normalization at a subsampling depth of 4000**
```bash
qiime diversity core-metrics-phylogenetic \
  --i-table filtered-table.qza \
  --i-phylogeny tree-no-chimera.qza \
  --p-sampling-depth 4000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results-4000
```
# Alpha diversity analysis
## Chao1
Chao1 alpha diversity analysis estimates species richness in a sample. It focuses on rare ASVs (appearing only once or twice), helping to assess how many species might be undetected due to limited sequencing depth. This is useful for comparing microbial richness across sample groups.
```bash
qiime diversity alpha \
  --i-table filtered-table.qza \
  --p-metric chao1 \
  --o-alpha-diversity chao1-diversity.qza
```
```bash
qiime metadata tabulate \
  --m-input-file chao1-diversity.qza \
  --o-visualization chao1-diversity.qzv
```
![image](https://github.com/user-attachments/assets/41a865d0-0e3e-4cf3-badb-3b6a511bfa0d)

```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity chao1-diversity.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization chao1-group-significance.qzv
```
![image](https://github.com/user-attachments/assets/420ce7d1-21bd-45b6-9258-c629552e0223)
The Kruskal-Wallis test results indicate:
    *Overall comparison (all groups)*
        H = 7.43, p-value = 0.115 → No significant difference in Chao1 alpha diversity among the groups (p > 0.05).
    *Pairwise comparisons*
        Most group comparisons have p-values > 0.05, suggesting no statistically significant differences.
        However, RP-20 vs RP-40 and RP-20 vs RP-5 have p-values = 0.0495, indicating a potential difference.
        Yet, the q-values are > 0.05 (after multiple testing correction), meaning the observed differences may not be strong enough to be considered statistically significant.
**Conclusion**: There is no significant difference in Chao1 alpha diversity among the groups after multiple comparison correction.

## Shannon index
The Shannon index measures alpha diversity, accounting for both species richness (number of species) and evenness (distribution of species abundances).
    Higher Shannon index → More diverse and evenly distributed microbial community.
    Lower Shannon index → A community dominated by a few species, indicating lower diversity.
It helps assess how microbial diversity changes across different conditions or treatment groups.j








