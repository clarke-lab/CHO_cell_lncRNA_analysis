#!/usr/bin/env bash
#### run transdecoder to predict the longest ORF from each feelnc transcript
#### inputs are: 1) reference genome sequence and 2) feelnc transcript gtf
#### 3) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-f = Reference genome sequence "
        echo "-g = feelnc gtf"
        echo "-o = Output directory"
        exit 2
fi
while getopts f:g:o: option
  do
    case "${option}"
      in
      f) REF_SEQ=${OPTARG};;
      g) FEELNC_GTF=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl \
$FEELNC_GTF $REF_SEQ > $OUT_DIR/candidate_lncrna_cdna.fa

~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs \
-t $OUT_DIR/candidate_lncrna_cdna.fa \
-S \
-O $OUT_DIR
