#!/usr/bin/env bash
#### StringTie transcriptome assembly
#### inputs are: 1) sample ID and 2) the inuput data directory 3) num threads
####  4) the reference GTF 5) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input bam directory"
        echo "-p = num processors"
        echo "-g = GTF file"
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
      g) GTF=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done

if [ ! -d $OUTDIR/individual_gtfs ]; then
mkdir -p $OUTDIR/individual_gtfs
fi

/mnt/HDD2/colin/bin/stringtie-2.0.3.Linux_x86_64/stringtie \
-p $THREADS \
-G $GTF \
--rf \
-j 5 \
-o $OUTDIR/individual_gtfs/"$SAMPLEID".gtf \
$INDIR/"$SAMPLEID"Aligned.sortedByCoord.out.bam

# END
