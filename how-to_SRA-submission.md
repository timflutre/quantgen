<!-- pandoc how-to_SRA-submission.md -f commonmark -t html -s -o how-to_SRA-submission.html -->

# How to submit to the NCBI's SRA?

*author: Timothée Flutre (INRA)*

*version: 06/09/2018*

Goal: insert sequencing reads to the [SRA](https://www.ncbi.nlm.nih.gov/sra) database of the [NCBI](https://www.ncbi.nlm.nih.gov).

All relevant, up-to-date information is in the official [doc](https://www.ncbi.nlm.nih.gov/sra/docs/submit/): this is only a quick tutorial!

## Prerequisites

### Raw data

List all your [fastq](https://en.wikipedia.org/wiki/FASTQ_format) files (compressed with [gzip](https://en.wikipedia.org/wiki/Gzip)) in a single file:

```
mkdir upload_NCBI
cd upload_NCBI
ls /my/data/*.fastq.gz > list_files_NCBI-SRA.txt
```

Number of files:
```
cat list_files_NCBI-SRA.txt | wc -l
```

Number of samples:
```
cat list_files_NCBI-SRA.txt | awk '{split($0, a, "_R"); print a[1]}' | sort | uniq -c | wc -l
```

### Meta-data

Make sure you have the passport data of each sample: unique identifier/code of the sample or the genotype, unique identifier/code of the cultivar, the sample names, species, etc.

For grapevine at INRA, we should use the [database of the RFCV](https://bioweb.supagro.inra.fr/collections_vigne/Home.php).

## Create an account

Go [here](https://www.ncbi.nlm.nih.gov/myncbi/).

## Create a BioProject submission

Go [here](https://submit.ncbi.nlm.nih.gov/subs/bioproject) and follow the steps.

Example:
* enter your Project title (e.g. "RAD-seq of ABC");
* enter your Public description;
* enter your Relevance (e.g. "agricultural");
* enter your Grants (e.g. "2014-ABC ; CASDAR");
* enter your Biosample type (e.g. "Plant");
* skip BioSamples for the moment;
* finish.

The status of your submission is available in your [submission portal](https://submit.ncbi.nlm.nih.gov/subs/).
Once processed, retrieve there your BioProject identifier, e.g PRJNA725341.

## Create a BioSample submission

Go [here](https://submit.ncbi.nlm.nih.gov/subs/biosample/) and follow the steps.

Advice: download the `.tsv` file and fill it programatically with your favorite programming language (R, Python, Perl, etc).

Notes:
* one can't submit a file with more than 1000 samples;
* the file should be in [ASCII](https://en.wikipedia.org/wiki/ASCII);
* one can't submit a sample already existing in the NCBI (the submission interface automatically detects such cases).

The status of your submission is available in your [submission portal](https://submit.ncbi.nlm.nih.gov/subs/).
Once processed, download the `attributes.tsv` file with BioSample accessions.

## Create a SRA submission

Go [here](https://submit.ncbi.nlm.nih.gov/subs/sra/) and follow the steps.

Download the `.tsv` file and fill it programatically with your favorite programming language (R, Python, Perl, etc).

Notes:
* make symbolic links of all fastq files into the same directory to ease the upload;
* use the Aspera Command-Line tool ([FASP](https://en.wikipedia.org/wiki/Fast_and_Secure_Protocol)) instead of [FTP](https://en.wikipedia.org/wiki/File_Transfer_Protocol) to upload many samples because it is much faster.

```
ls fastq_gbs-279/* | wc -l
ls fastq_gbs-279/*_R*.fastq.gz | wc -l
ls -la fastq_gbs-279/ | grep "\->" | wc -l
nohup ascp -i ~/.ssh/aspera.openssh -QT -l100m -k1 \
      -d fastq_gbs-279 \
      subasp@upload.ncbi.nlm.nih.gov:uploads/your@email_<...> \
      >& stdout_ascp.txt &
```
