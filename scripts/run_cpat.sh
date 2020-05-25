#!/usr/bin/env bash
#### determine the protein coding potential of putative lncRNA transcripts
#### inputs are: 1) fasta file of transcripts and 2) the output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-f = lncRNA fasta sequence"
        echo "-o = Output directory"
        exit 2
fi
while getopts f:o: option
  do
    case "${option}"
      in
      f) LNCRNA_SEQ=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

#get mouse CPAT data if needed
if [ ! -f $OUT_DIR/Mouse_logitModel.RData ]; then
wget https://ayera.dl.sourceforge.net/project/rna-cpat/v1.2.2/prebuilt_model/Mouse_logitModel.RData \
-P $OUT_DIR
fi

if [ ! -f $OUT_DIR/Mouse_Hexamer.tsv ]; then
wget https://ayera.dl.sourceforge.net/project/rna-cpat/v1.2.2/prebuilt_model/Mouse_Hexamer.tsv \
-P $OUT_DIR
fi

# run CPAT
cpat.py \
-g $LNCRNA_SEQ \
-o $OUT_DIR/CPAT.analysis.txt \
-x $OUT_DIR/Mouse_Hexamer.tsv \
-d $OUT_DIR/Mouse_logitModel.RData

# END
