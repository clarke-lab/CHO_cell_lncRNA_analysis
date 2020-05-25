#!/usr/bin/env bash
if (($# == 0)); then
        echo "Usage:"
        echo "-G = Reference genome GTF"
        echo "-g = Stringtie GTF"
        echo "-o = Output directory"
        exit 2
fi
while getopts f:G:g:o: option
  do
    case "${option}"
      in
      G) REFGTF=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done

# CREATE A NEW GTF USING NOVEL AND PREVIOUSLY ANNOTATED LNCRNAS FROM FEELnc_codpot GTF
grep -wFf $OUTDIR/novel_lncrna_transcripts.list $OUTDIR/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf > $OUTDIR/novel_lncrna.gtf
grep -wFf reference_genome/ensembl_lncrna_transcript.list $OUTDIR/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf > $OUTDIR/ensembl_lncrna.gtf
cat $OUTDIR/novel_lncrna.gtf $OUTDIR/ensembl_lncrna.gtf > $OUTDIR/all_lncrna.gtf


awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $REFGTF | grep 'protein_coding' > $OUTDIR/known_mrna.gtf

#FEELnc_classifier.pl \
#-i $OUTDIR/all_lncrna.gtf \
#-a $OUTDIR/known_mrna.gtf \
#> $OUTDIR/lncRNA_classes.txt

Rscript R/filter_lncrna.R "lncrna_annotation/lncRNA_classes.txt" "lncrna_annotation"
