#!/usr/bin/env bash
## anotation of lncRNAs via RFAM blast and syntheny with human and mouse

if (($# == 0)); then
        echo "Usage:"
        echo "-s lncRNA sequences"
        echo "-o = Output directory"
        exit 2
fi
while getopts g:o:l:k:a: option
  do
    case "${option}"
      in
      a) STR_GTF=${OPTARG};;
      k) ENS_LNCRNA=${OPTARG};;
      l) LNCRNA_GTF=${OPTARG};;
      g) REF_GTF=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

# each lncrna is classified using FEELNc
mkdir $OUT_DIR/classification
FEELnc_classifier.pl   -i $LNCRNA_GTF -a $REF_GTF > $OUT_DIR/classification/first_classification.txt

#awk '{ if ($3 == "exon") print $0; }' lncrna_annotation/all_lncrna.gtf | grep -o "transcript_id [^;]\+;" | awk -v FS=" " '{ arr[substr($2, 2, length($2)-3)] += 1; } END { for (t in arr) { if (arr[t] > 1) { print t;} } }'

# all monexonic
mkdir $OUT_DIR/monoexonic_filter

# determine which of the putative lncRNAs have 1 exon
sed 's/\s/\t/g'  $LNCRNA_GTF | awk '$3 == "exon"'| datamash -s  -g 12 count 3 \
| awk '$2 < 2' |  tr -d \" | tr -d \; | awk '{print $1}' | sort > $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt

sed 's/\s/\t/g' $LNCRNA_GTF | awk '$3 == "exon"'| datamash -s  -g 12 count 3 \
| awk '$2 >= 2' |   tr -d \" | tr -d \; | awk '{print $1}' | sort > $OUT_DIR/monoexonic_filter/multiexonic_lncrnas.txt

# one exon lncRNAs are overlapped with the previous RFAM blast
# determine monoexonic with hits in rfam, as well as synteny with human mouse
# rfam blast hits
awk '{print $1}' $OUT_DIR/RFAM/lncrna.blastn.outfmt6 | sort > \
$OUT_DIR/monoexonic_filter/rfam_hit_transcripts.txt

comm -12 $OUT_DIR/monoexonic_filter/rfam_hit_transcripts.txt $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt >  \
$OUT_DIR/monoexonic_filter/rfam_monoexonic_lncrnas.txt


awk '{print $1}' $OUT_DIR/liftover/lncrna_cho_to_human_otho.txt | sort > $OUT_DIR/monoexonic_filter/human_synteny.txt
comm -12 $OUT_DIR/monoexonic_filter/human_synteny.txt $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt >  \
$OUT_DIR/monoexonic_filter/human_monoexonic_lncrnas.txt

awk '{print $1}' $OUT_DIR/liftover/lncrna_cho_to_mouse_ortho.txt | sort  > lncrna_annotation/monoexonic_filter/mouse_synteny.txt
comm -12 $OUT_DIR/monoexonic_filter/mouse_synteny.txt $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt >  \
$OUT_DIR/monoexonic_filter/mouse_monoexonic_lncrnas.txt

# retain monoexonic classified as antisense to a protein coding gene
awk '$8 ==0 {print}' lncrna_annotation/classification/first_classification.txt | grep 'antisense' | awk '{print $3}' | uniq | sort > $OUT_DIR/monoexonic_filter/feelnc_antisense_overlap.txt
comm -12 $OUT_DIR/monoexonic_filter/feelnc_antisense_overlap.txt $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt > \
$OUT_DIR/monoexonic_filter/antisense_monoexonic_lncrnas.txt

# ensembl annotated lncrnas

# keep lnrnas 1) monoexonic found by blast against rfam 2) synteny,
# 3) antisense to protein coding genes 4) annotated in ensembl and with more than 1 exon
cat \
reference_genome/ensembl_lncrna_transcript.txt \
$OUT_DIR/monoexonic_filter/multiexonic_lncrnas.txt \
$OUT_DIR/monoexonic_filter/rfam_monoexonic_lncrnas.txt \
$OUT_DIR/monoexonic_filter/human_monoexonic_lncrnas.txt \
$OUT_DIR/monoexonic_filter/mouse_monoexonic_lncrnas.txt \
$OUT_DIR/monoexonic_filter/antisense_monoexonic_lncrnas.txt \
$ENS_LNCRNA |  sort| uniq  > lncrna_annotation/monoexonic_filter/final_lncrna_list.txt

grep -wFf $OUT_DIR/monoexonic_filter/final_lncrna_list.txt $STR_GTF > $OUT_DIR/lncrna.gtf


## rerun classifer on final lncrna list
FEELnc_classifier.pl \
-i $OUT_DIR/lncrna.gtf \
-a $REF_GTF \
> $OUT_DIR/classification/second_classification.txt
