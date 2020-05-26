#!/usr/bin/env bash
#### calculate the TPM expression for each transcript
#### in the stringtie assembly for each sample
#### inputs are: 1) sample and 2) stringtie gtf 3) bam file location
#### 4) processors 5) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = Sample ID "
        echo "-g = Stringtie GTF"
        echo "-b = star bam directory"
        echo "-p = num processors"
        echo "-o = Output directory"
        exit 2
fi
while getopts s:b:g:p:o: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      b) BAM_DIR=${OPTARG};;
      g) STR_GTF=${OPTARG};;
      p) THREADS=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR/"$SAMPLEID"
fi


/mnt/HDD2/colin/bin/stringtie-2.0.3.Linux_x86_64/stringtie \
-p $THREADS \
-G $STR_GTF \
-e \
-B \
-o $OUT_DIR/"$SAMPLE_ID"/transcripts.gtf \
-A $OUT_DIR/"$SAMPLE_ID"/gene_abundances.tsv \
$BAM_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam
