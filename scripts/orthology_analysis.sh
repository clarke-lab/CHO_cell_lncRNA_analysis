#!/usr/bin/env bash
#### anotation of lncRNAs via RFAM and syntheny with human and mouse
#### inputs are: 1-3)output directories 4) evalue cutoff 5) transcript sequences
#### 6) protein sequences 6) threads
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

#!/usr/bin/env bash
## anotation of lncRNAs via RFAM blast and syntheny with human and mouse

if (($# == 0)); then
        echo "Usage:"
        echo "-s lncRNA sequences"
        echo "-o = Output directory"
        exit 2
fi
while getopts s:g:o: option
  do
    case "${option}"
      in
      s) lncRNA_fasta=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done

if [ ! -d $OUTDIR/liftover ]; then
mkdir -p $OUTDIR/liftover
fi

# download UCSC version of the chok1 genome and chain files for CHOK1
if [ ! -f $OUTDIR/liftover/criGriChoV1.fa ] || [ ! -f $OUTDIR/liftover/criGriChoV1ToMm10.over.chain ] || [ ! -f $OUTDIR/liftover/criGriChoV1ToHg38.over.chain.gz ] ; then
  wget http://hgdownload.soe.ucsc.edu/goldenPath/criGriChoV1/bigZips/criGriChoV1.fa.gz \
  -P $OUTDIR/liftover/
  wget http://hgdownload.soe.ucsc.edu/goldenPath/criGriChoV1/liftOver/criGriChoV1ToHg38.over.chain.gz -P $OUTDIR/liftover/
  wget http://hgdownload.soe.ucsc.edu/goldenPath/criGriChoV1/liftOver/criGriChoV1ToMm10.over.chain.gz -P $OUTDIR/liftover/
  gunzip $OUTDIR/liftover/*.gz
fi

#build gmap index for the genome
if [ ! -d $OUTDIR/liftover/UCSC ]; then
echo "building gmap index - can take some time"
#gmap_build -D $OUTDIR/liftover/ -d UCSC $OUTDIR/liftover/criGriChoV1.fa
fi

#map lncRNA sequences to UCSC genome and output gff3 file
gmap -D $OUTDIR/liftover/UCSC \
     -d UCSC \
     -n 1 \
     --no-chimeras \
     -f gff3_gene \
     -t 32 \
     $lncRNA_fasta > $OUTDIR/liftover/gmap.gff3

#convert the gff3 from gmap to bed
awk -F'[\t;]' 'OFS="\t" {if ($3=="gene"){print $1,$4-1,$5,$10,".",$7}}' $OUTDIR/liftover/gmap.gff3 | sed 's/Name=//g' >  \
$OUTDIR/liftover/gmap.bed

# use the chain files to do the liftover for human and mouse
/mnt/HDD2/colin/bin/liftOver -gff -minMatch=0.1 \
$OUTDIR/liftover/gmap.gff3 \
$OUTDIR/liftover/criGriChoV1ToHg38.over.chain \
$OUTDIR/liftover/choTohuman.bed \
$OUTDIR/liftover/unmap.human.bed

/mnt/HDD2/colin/bin/liftOver -gff -minMatch=0.1 \
$OUTDIR/liftover/gmap.gff3 \
$OUTDIR/liftover/criGriChoV1ToMm10.over.chain \
$OUTDIR/liftover/choTomouse.bed \
$OUTDIR/liftover/unmap.mouse.bed

# get the gencode GTF files and covert to bed
if [ ! -f $OUTDIR/liftover/gencode.vM23.long_noncoding_RNAs.gtf ] || [ ! -f $OUTDIR/liftover/gencode.v32.long_noncoding_RNAs.gtf ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.long_noncoding_RNAs.gtf.gz \
  -P $OUTDIR/liftover
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz \
  -P $OUTDIR/liftover
  gunzip $OUTDIR/liftover/*.long_noncoding_RNAs.gtf.gz
fi

# convert gencode gtfs to bed and intersect
awk 'OFS="\t" {if ($3=="gene"){print $1,$4-1,$5,$10,".",$7}}' $OUTDIR/liftover/gencode.v32.long_noncoding_RNAs.gtf|\
sed 's/;//g;s/"//g' >  $OUTDIR/liftover/gencode.v32.long_noncoding_RNAs.bed
bedtools intersect -wb -a "$OUTDIR"/liftover/choTohuman.bed -b $OUTDIR/liftover/gencode.v32.long_noncoding_RNAs.bed > \
$OUTDIR/liftover/human.conserved.bed

awk 'OFS="\t" {if ($3=="gene"){print $1,$4-1,$5,$10,".",$7}}' $OUTDIR/liftover/gencode.vM23.long_noncoding_RNAs.gtf|\
sed 's/;//g;s/"//g' >  $OUTDIR/liftover/gencode.vM23.long_noncoding_RNAs.bed
bedtools intersect -wb -a "$OUTDIR"/liftover/choTomouse.bed -b $OUTDIR/liftover/gencode.vM23.long_noncoding_RNAs.bed > \
$OUTDIR/liftover/mouse.conserved.bed

#synteny annotation of CHOK1 lncRNAs
grep -wFf $OUTDIR/firstpass_filter/all_lncrna_transcripts.txt lncrna_annotation/liftover/human.conserved.bed \
| grep 'gene' | awk '{print $9,$13}' | awk -F= '{print $3}' | sed -e 's/\.[0-9]//2' | sed -e 's/\.[0-9]//2'   > $OUTDIR/liftover/lncrna_cho_to_human.txt.tmp
awk '{gsub(/\"|\.*/,"",$2)}1' $OUTDIR/liftover/lncrna_cho_to_human.txt.tmp  > $OUTDIR/liftover/lncrna_cho_to_human.txt

grep -wFf $OUTDIR/firstpass_filter/all_lncrna_transcripts.txt lncrna_annotation/liftover/mouse.conserved.bed \
| grep 'gene' | awk '{print $9,$13}' | awk -F= '{print $3}' | sed -e 's/\.//2' | sed -e 's/\.//2' > $OUTDIR/liftover/lncrna_cho_to_mouse.txt.tmp
awk '{gsub(/\"|\.*/,"",$2)}1' $OUTDIR/liftover/lncrna_cho_to_mouse.txt.tmp > $OUTDIR/liftover/lncrna_cho_to_mouse.txt

grep -wFf $OUTDIR/firstpass_filter/all_lncrna_transcripts.txt lncrna_annotation/liftover/human.conserved.bed \
| grep 'gene' | awk '{print $9,$13}' | awk -F= '{print $3}'  > $OUTDIR/liftover/lncrna_cho_to_human_otho.txt

grep -wFf $OUTDIR/firstpass_filter/all_lncrna_transcripts.txt lncrna_annotation/liftover/mouse.conserved.bed \
| grep 'gene' | awk '{print $9,$13}' | awk -F= '{print $3}'  > $OUTDIR/liftover/lncrna_cho_to_mouse_ortho.txt
