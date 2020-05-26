#!/usr/bin/env bash
#### count reads aligned to reference genome for differential expression
#### inputs are: 1) sample ID and 2) bam input directory 2) reference GTF
#### 4) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019
if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-m = Mapping directory"
        echo "-g = GTF"
        echo "-o = path for output BAMs"
        exit 2
fi
while getopts s:m:g:o: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      g) REF_GTF=${OPTARG};;
      m) MAP_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

mkdir -p $OUT_DIR
echo $SAMPLE_ID
htseq-count -r pos -f bam -i gene_id -s reverse $MAP_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam $REF_GTF > "$OUT_DIR"/"$SAMPLE_ID".counts

# END
