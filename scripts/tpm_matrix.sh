#!/usr/bin/env bash
#### merge the TPM values for each sample into a single matix
#### inputs are: 1) a list of sample IDs and 2) a storage directory for
#### for the results
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID "
        echo "-o = Output directory"
        exit 2
fi
while getopts s:o: option
  do
    case "${option}"
      in
      s) SAMPLE_LIST=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

input_directories=$(sed "s|^|$OUT_DIR|" $SAMPLE_LIST | paste -sd ",")


./scripts/make_tpm_expression_matrix.pl \
--expression_metric=TPM \
--result_dirs=$input_directories \
--transcript_matrix_file=$OUT_DIR/transcript_tpm_all_samples_v1.tsv \
--gene_matrix_file=$OUT_DIR/gene_tpm_all_samples.tsv

# clean "na" from the matrices
grep -v na $OUT_DIR/transcript_tpm_all_samples_v1.tsv > $OUT_DIR/transcript_tpm_all_samples.tsv
rm $OUT_DIR/transcript_tpm_all_samples_v1.tsv
