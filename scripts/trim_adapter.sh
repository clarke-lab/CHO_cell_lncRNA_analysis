#!/usr/bin/env bash
#### Trim adapter sequences using cutadapt
#### inputs are: 1) sample ID and 2) the data directory & output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-o = output directory"
        exit 2
fi

while getopts s:i:o: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      i) IN_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

# trim adapter using cutadapt
# Following https://cutadapt.readthedocs.io/en/stable/guide.html
# Using the longer adapter sequences recommended by this document
cutadapt  \
-A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
-o $OUT_DIR/"$SAMPLE_ID"_1.fastq.gz \
-p $OUT_DIR/"$SAMPLE_ID"_2.fastq.gz \
   $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz

# END
