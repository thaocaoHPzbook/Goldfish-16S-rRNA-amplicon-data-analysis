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
# Taxonomic clasification
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



