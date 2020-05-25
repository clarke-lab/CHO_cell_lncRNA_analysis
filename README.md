# CHO cell alternative splicing analysis

## Introduction
---
This respository enables the reproduction of the analysis described in:

Motheramgari *et al.* **Expanding the Chinese hamster ovary cell long non-coding RNA transcriptome using RNASeq**
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
In addition, a complementary annotation file is created from NCBI to help with annotation later
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

## lncRNA annotation
Make a directory to hold each of the different analyses
```bash
mkdir lncrna_annotation
```
### Calculate the expression values for each transcript
```bash
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
  ./scripts/calc_tpm.sh -p 32 -s $sample -g transcriptome_assembly/stringtie.all.transcripts.gtf -o lncrna_annotation/TPM -b data/mapped
done
```
```bash
# list samples
awk 'NR>1 {print $2}' data/sample_info.txt > lncrna_annotation/TPM/sample_list.txt
./scripts/tpm_matrix.sh -s lncrna_annotation/TPM/sample_list.txt -o lncrna_annotation/TPM/
```

### FEELNc analysis
The stringtie transcriptome assembly is used to predict lncRNAs using FEELNc
```bash
./scripts/run_feelnc.sh -G reference_genome/ensembl_chok1_genome.gtf -g transcriptome_assembly/non_protein_coding_stringtie.gtf -f reference_genome/ensembl_chok1_genome.fa -o lncrna_annotation/FEELnc
```

### Assess FEELNc output using additional protein potential calculators, PFAM search and BLAST against protein and RNA databases
#### Use transdecoder to create a cDNA fasta file identify the longest ORF for each candidate lncRNA
```bash
./scripts/run_transdecoder.sh -g lncrna_annotation/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -f reference_genome/ensembl_chok1_genome.fa -o lncrna_annotation/TRANSDECODER
```
#### CPAT coding prediction for FEELNc candidate lncRNAs
```bash
./scripts/run_cpat.sh -f lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa -o lncrna_annotation/CPAT
```
#### CPC2 coding prediction for FEELnc candidate lncRNAs
```bash
  ./scripts/run_cpc2.sh -f lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa -o lncrna_annotation/CPC2
```
#### Assess FEELnc candiate lncRNAs for the presence of protein domains
```bash
./scripts/run_hmmscan.sh -t 32 -e 1e-5 -p lncrna_annotation/TRANSDECODER/longest_orfs.pep -o lncrna_annotation/PFAM
```
#### Assess FEELnc candiate lncRNAs for the presence of proteins, miRNAs, and other non-coding RNAs (e.g. snoRNAs) using BLAST
```bash
./scripts/run_blast.sh  -t 32 -e 1e-5 -s lncrna_annotation/SWISSPROT -p lncrna_annotation/TRANSDECODER/longest_orfs.pep -m lncrna_annotation/MIRBASE -r lncrna_annotation/RFAM -n lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa
```

### intersect the results
```bash
/usr/bin/Rscript R/filter_lncrna.R \
  "lncrna_annotation/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf" \
  "lncrna_annotation/CPC2/CPC2.analysis.txt" \
  "lncrna_annotation/CPAT/CPAT.analysis.txt" \
  "lncrna_annotation/SWISSPROT/blastp.outfmt6" \
  "lncrna_annotation/MIRBASE/blastn.outfmt6" \
  "lncrna_annotation/PFAM/pfam_domain_transcripts" \
  "lncrna_annotation/RFAM/blastn.outfmt6" \
  "lncrna_annotation/TPM/transcript_tpm_all_samples.tsv" \
  "lncrna_annotation/FEELnc/lncRNA_classes.txt" \
  "transcriptome_assembly/non_protein_coding_stringtie.gtf" \
  "lncrna_annotation"
```

### Second stage of filtering
Here we filter monoexonic annotations. First we determine the lncRNAs that pass the first stage of filtering and their
orthology with human and mouse lncRNAs
#### determine synteny with human and mouse gencode lncRNAs
```bash
./scripts/orthology_analysis.sh -s lncrna_annotation/firstpass_filter/lncRNA.fasta -o lncrna_annotation
```
## Filter monoexonic lncRNAs

#### keep only monoexonic lncRNAs that are:
#### 1) antisense to a protein coding gene 2) are orthologous with human or mouse 3) annotated as lncRNA in ensembl
```bash
./scripts/filter_monoexonic_lncrnas.sh \
-o lncrna_annotation \
-g reference_genome/ensembl_chok1_protein.gtf \
-l lncrna_annotation/firstpass_filter/all_lncrna.gtf \
-k reference_genome/ensembl_lncrna_transcript.txt \
-a transcriptome_assembly/non_protein_coding_stringtie.gtf
```
#
```bash
/usr/bin/Rscript R/Rscript R/simplify_class.R \
  "lncrna_annotation/classification/second_classification.txt" \
  "lncrna_annotation/classification/final_classification.txt"
```

## Differential expression analysis
#### Gene level counting
```bash
mkdir differential_expression
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/htseq_count.sh -s $sample -m data/mapped -g transcriptome_assembly/stringtie_original.appended.fp.s.filtered.gtf -o differential_expression/counts&
done
```

```bash
mkdir differential_expression/results
Rscript R/run_deseq.R "differential_expression/counts" "differential_expression/results"
```
