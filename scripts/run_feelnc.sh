#!/usr/bin/env bash
#### Carry out initial lncRNA prediction using feelnc
#### inputs are: 1) reference genome fasta 2) reference genome GTF
#### 3) stringtie GTF 4) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-f = Reference genome sequence "
        echo "-G = Reference genome GTF"
        echo "-g = Stringtie GTF"
        echo "-o = Output directory"
        exit 2
fi
while getopts f:G:g:o: option
  do
    case "${option}"
      in
      f) REF_SEQ=${OPTARG};;
      G) REF_GTF=${OPTARG};;
      g) STR_GTF=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

# filter transcripts overlapping with sense protein coding exons.
# Keeping monoex to deal with misc_RNA biotype
FEELnc_filter.pl \
-i $STR_GTF \
-a $REF_GTF \
-b transcript_biotype=protein_coding,pseudogene  \
--monoex=1 \
--size=200 \
-p 32 \
> $OUT_DIR/candidate_lncRNA.gtf

# create a gtf for known protein coding transcripts
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' \
$REF_GTF | \
grep 'protein_coding' \
> $OUT_DIR/known_mrna.gtf

#determine protein coding potential
FEELnc_codpot.pl \
-i $OUT_DIR/candidate_lncRNA.gtf \
-m shuffle \
-a $OUT_DIR/known_mrna.gtf \
-g $REF_SEQ \
--verbosity=0 \
--outdir $OUT_DIR/feelnc_codpot_out/
