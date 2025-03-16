# 1. Dowload data and create the list of sample's name
## Create the Goldfish directory with a subdirectory named Raw_data
```bash
mkdir -p Home/hp/Goldfish/Raw_data
```
## Navigate to the Raw_data directory
```bash
cd Home/hp/Goldfish/Raw_data
```

## Download Goldfish Rawdata file into the Raw_data directory
```bash
wget -P Home/hp/Goldfish/Raw_data <URL_FILE>
```

## Find all .fastq files in Goldfish/Raw_data and save their filenames to IDs.list
```bash
find Home/hpGoldfish/Raw_data -name "*.fastq" | sed 's|.*/||' > Home/hp/Goldfish/Raw_data/IDs.list
```
## Set read permissions for the IDs.list file
```bash
chmod +r Home/hpGoldfish/Raw_data/IDs.list
```
# 2. Quality check
## Create script file to run fastqc
```bash
nano run_fastqc.sh
```
The script file named [run_fastqc.sh](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Quality_Check/run_fastqc.sh)
 is generated.

```bash
chmod +x run_fastqc.sh
./run_fastqc.sh
```
QC results will be genenerate in folder named **fastqc_results**, find **[multiqc_report.html](https://github.com/thaocaoHPzbook/Goldfish-16S-rRNA-amplicon-data-analysis/blob/main/Quality_Check/multiqc_report.html)** file to read the results.

![image](https://github.com/user-attachments/assets/01848d8e-5dfc-4298-baa5-73a32b2490fb)

The MultiQC results indicate that the data quality is quite good, with all quality metrics meeting the expected standards. Additionally, the adapter contents have been completely removed, ensuring clean and reliable sequences for downstream analysis.



