# CHO cell alternative splicing analysis

## Introduction
---
This respository enables the reproduction of the analysis described in:

Tzani *et al.* **Sub physiological temperature induces pervasive alternative splicing in Chinese hamster ovary cells**
*bioRxiv 2019* [https://doi.org/10.1101/863175](https://doi.org/10.1101/863175)

### Data availability:
[NIBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA593052/)  
[European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJNA593052)

### Dependancies
All the programmes must be added to the PATH to run the workflow
- [Python 2.7.12](https://www.python.org/download/releases/2.7/)
- [trimmomatic 0.36](http://www.usadellab.org/cms/?page=trimmomatic)
- [cutadpat 1.18](https://cutadapt.readthedocs.io/en/stable/)
- [STAR-2.7.2d](https://github.com/alexdobin/STAR)
- [rMats.4.0.2](http://rnaseq-mats.sourceforge.net)
- [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot)
- [stringtie 2.0.3](http://ccb.jhu.edu/software/stringtie/index.shtml)
- [samtools 1.6](http://www.htslib.org)

- R 3.5.2
    - [dpylr 0.8.4](https://dplyr.tidyverse.org)
    - [DESeq2 1.22.0](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    - [xlxs 0.6.2](https://cran.r-project.org/web/packages/xlsx/index.html)
    - [stringr 1.4.0](https://stringr.tidyverse.org)
    - [biomaRt 2.38.0](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

## RNASeq data preprocesssing
---
### Download the data from ENA
This is a simple way to dowload from ENA, for higher speed download use the Aspera client
Total data download size: **~95G**
```bash
mkdir -p data/ena
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/057/SRR10572657/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/058/SRR10572658/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/059/SRR10572659/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/060/SRR10572660/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/061/SRR10572661/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/062/SRR10572662/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/063/SRR10572663/*" -P data/ena || { handle ; error ; }

wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/064/SRR10572664/*" -P data/ena || { handle ; error ; }
```
### trim adapter sequences
Adapter trimming for the Illumina TruSeq adapters
```bash
mkdir data/cutadapt
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/trim_adapter.sh -s $sample  -i data/ena -o data/cutadapt
done
```
### quality trimming
Trimmomatic for removing low quality bases and filtering resulting reads that are  
too short. The script makes two subfolders in the trimmomatic folder for paired and unpaired reads  
```bash
mkdir data/preprocessed
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/trim_quality.sh -s $sample -i data/cutadapt -o data/preprocessed -p 32
done
```
## Read Mapping
### Download the reference genome and NCBI annotation for CHO K1
. The sequence is also prepped for further analysis by retaining only the scaffold ID.
In addition, a complementay annotation file is created from NCBI to help with annotation later
```bash
mkdir reference_genome
./scripts/prepare_genome.sh -v 98 -o reference_genome
```
### make the STAR index
Use the reference genome to build an index for mapping the RNASeq reads
```bash
./scripts/make_star_index.sh -g reference_genome/ensembl_chok1_genome.fa -a reference_genome/ensembl_chok1_genome.gtf -p 32 -d reference_genome
```

### map to CHOK1 genome
```bash
mkdir data/mapped
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/star_mapping.sh -s $sample -i data/preprocessed/paired -g reference_genome/star_index -o data/mapped -p 32
done
```

## Genome guided assembly
### string tie assembly
```bash
mkdir transcriptome_assembly
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
  ./scripts/run_stringtie.sh -s $sample -i data/mapped -g reference_genome/ensembl_chok1_genome.gtf -o transcriptome_assembly -p 32
done
```

### merge individual stringtie assemblies and compare to ENSEMBL annotation
```bash
./scripts/stringtie_merge.sh -t transcriptome_assembly -g reference_genome/ensembl_chok1_genome.gtf -r reference_genome
```
