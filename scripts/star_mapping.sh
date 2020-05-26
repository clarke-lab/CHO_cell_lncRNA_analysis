#!/usr/bin/env bash
#### Map reads to the reference index
#### inputs are: 1) sample ID and 2) the data directory 3) output directory
####  4) number of threads
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-p = num processors"
        echo "-g = path to star genome index"
        echo "-o = path for output BAMs"
        exit 2
fi
while getopts s:i:p:g:o: option
  do
    case "${option}"
      in
      s) SAMPLEID=${OPTARG};;
      i) INDIR=${OPTARG};;
      p) THREADS=${OPTARG};;
      g) INDEX=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done


/mnt/HDD2/colin/bin/STAR-2.7.2d/bin/Linux_x86_64/STAR \
--runThreadN $THREADS \
--readFilesIn $INDIR/"$SAMPLEID"_1.fastq.gz $INDIR/"$SAMPLEID"_2.fastq.gz \
--genomeDir $INDEX \
--readFilesCommand gunzip -c \
--outFileNamePrefix $OUTDIR/"$SAMPLEID" \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic

#--outSAMstrandField intronMotif \

# create a BAM index
samtools index $OUTDIR/"$SAMPLEID"Aligned.sortedByCoord.out.bam

# END
