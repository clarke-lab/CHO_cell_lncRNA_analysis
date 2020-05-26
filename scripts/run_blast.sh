#!/usr/bin/env bash
#### blast against RFAM,SWISSPROT and MIRBASE with feelnc transcripts and
#### transdecoder sequences
#### inputs are: 1-3)output directories 4) evalue cutoff 5) transcript sequences
#### 6) protein sequences 6) threads
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = SWISSPROT BLAST output dir"
        echo "-r = RFAM BLAST output dir"
        echo "-m = MIRBASE BLAST output dir"
        echo "-e = E-Value threshold"
        echo "-n = FEELNc lncRNA nucleotide sequence"
        echo "-p = protein sequence"
        echo "-t = Processor numbers"
        exit 2
fi
while getopts r:s:e:p:n:s:m:t: option
  do
    case "${option}"
      in
      r) RFAM_OUTPUT_DIR=${OPTARG};;
      e) EVALUE=${OPTARG};;
      p) PROTEIN=${OPTARG};;
      n) NUCLEOTIDE=${OPTARG};;
      s) SWISSPROT_OUTDIR=${OPTARG};;
      m) MIRBASEOUTDIR=${OPTARG};;
      t) THREADS=${OPTARG};;
    esac
done

# Swiss prot
if [ ! -d $SWISSPROT_OUTDIR ]; then
mkdir -p $SWISSPROT_OUTDIR
fi

if [ ! -f $SWISSPROT_OUTDIR/uniprot_sprot.fasta ]; then
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
-P $SWISSPROT_OUTDIR
gunzip $SWISSPROT_OUTDIR/uniprot_sprot.fasta.gz
fi

if [ ! -f $SWISSPROT_OUTDIR/uniprot_sprot.fasta.phr ]; then
makeblastdb -in $SWISSPROT_OUTDIR/uniprot_sprot.fasta -dbtype prot
fi

blastp \
-query $PROTEIN \
-db $SWISSPROT_OUTDIR/uniprot_sprot.fasta  \
-max_target_seqs 1 \
-outfmt 6 \
-evalue $EVALUE \
-num_threads $THREADS \
> $SWISSPROT_OUTDIR/blastp.outfmt6

#miRBase
if [ ! -d $MIRBASEOUTDIR ]; then
mkdir -p $MIRBASEOUTDIR
fi

if [ ! -f $MIRBASEOUTDIR/hairpin.fa ]; then
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz \
-P $MIRBASEOUTDIR
gunzip $MIRBASEOUTDIR/hairpin.fa.gz
fi

if [ ! -f $MIRBASEOUTDIR/hairpin.fa.nhr ]; then
makeblastdb -in $MIRBASEOUTDIR/hairpin.fa -dbtype nucl
fi

blastn \
-db $MIRBASEOUTDIR/hairpin.fa \
-query $NUCLEOTIDE \
-strand plus \
-evalue $EVALUE \
-num_threads $THREADS \
-outfmt 6 \
-num_alignments 1 \
> $MIRBASEOUTDIR/blastn.outfmt6

# RFAM
if [ ! -d $RFAM_OUTPUT_DIR ]; then
mkdir -p $RFAM_OUTPUT_DIR
fi

# create 2 files for blast snorna fasta for filtering here and lncRNA for comparison later
if [ ! -f $RFAM_OUTPUT_DIR/rfam.filter.fasta ] || [ ! -f $RFAM_OUTPUT_DIR/rfam_lncrna.fasta ] ; then
  wget -r ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/ -P $RFAM_OUTPUT_DIR
  awk {'print "lncrna_annotation/RFAM/ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$1".fa.gz"'} data/rfam_lncrna.info > $RFAM_OUTPUT_DIR/rfam_lncrna.fasta.txt
  awk {'print "lncrna_annotation/RFAM/ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$1".fa.gz"'} data/rfam_snorna.info > $RFAM_OUTPUT_DIR/rfam_snorna.fasta.txt
  awk {'print "lncrna_annotation/RFAM/ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$1".fa.gz"'} data/rfam_ribozyme.info > $RFAM_OUTPUT_DIR/rfam_ribozyme.fasta.txt

  #lncrna_fasta
  cat $RFAM_OUTPUT_DIR/rfam_lncrna.fasta.txt | while read rfam_fasta; do
      zcat $rfam_fasta >> $RFAM_OUTPUT_DIR/rfam_lncrna.fasta
    done

  cat $RFAM_OUTPUT_DIR/rfam_snorna.fasta.txt | while read rfam_fasta; do
     zcat $rfam_fasta >> $RFAM_OUTPUT_DIR/rfam_snorna.fasta
   done

   cat $RFAM_OUTPUT_DIR/rfam_ribozyme.fasta.txt | while read rfam_fasta; do
      zcat $rfam_fasta >> $RFAM_OUTPUT_DIR/rfam_ribozyme.fasta
    done

    rm -r $RFAM_OUTPUT_DIR/ftp.ebi.ac.uk
fi

# merge the RFAM sequences to filter against
cat $RFAM_OUTPUT_DIR/rfam_snorna.fasta $RFAM_OUTPUT_DIR/rfam_ribozyme.fasta > $RFAM_OUTPUT_DIR/rfam.filter.fasta
#create blastdb for RFAM lncRNA and non-coding RNA
if [ ! -f $RFAM_OUTPUT_DIR/rfam.filter.fasta.hnr ]; then
makeblastdb -in $RFAM_OUTPUT_DIR/rfam.filter.fasta -dbtype nucl
fi

if [ ! -f $RFAM_OUTPUT_DIR/rfam_lncrna.fasta.hnr ]; then
makeblastdb -in $RFAM_OUTPUT_DIR/rfam_lncrna.fasta -dbtype nucl
fi

blastn \
-db $RFAM_OUTPUT_DIR/rfam.filter.fasta \
-query $NUCLEOTIDE \
-strand plus \
-evalue $EVALUE \
-perc_identity 95 \
-num_threads $THREADS \
-outfmt 6 \
-num_alignments 1 \
> $RFAM_OUTPUT_DIR/blastn.outfmt6


blastn \
-db $RFAM_OUTPUT_DIR/rfam_lncrna.fasta \
-query $NUCLEOTIDE \
-strand plus \
-evalue $EVALUE \
-perc_identity 90 \
-num_threads $THREADS \
-outfmt 6 \
-num_alignments 1 \
> $RFAM_OUTPUT_DIR/lncrna.blastn.outfmt6

echo "----------------------"
echo "BLAST steps complete"
echo "----------------------"
echo "protein BLAST output written to $SWISSPROT_OUTDIR/blastp.outfmt6"
echo "mriRNA BLAST output written to $MIRBASEOUTDIR/blastn.outfmt6"
echo "snoRNA and ribozyme BLAST output written to $RFAM_OUTPUT_DIR/blastn.outfmt6"
echo "lncRNA BLAST output written to $RFAM_OUTPUT_DIR/lncrna.blastn.outfmt6"
