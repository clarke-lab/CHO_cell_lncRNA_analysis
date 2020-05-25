#!/usr/bin/env bash
#### make a reference index for STAR
#### inputs are: 1) ENSEMBL GTF and 2) ENSEMBL FASTA 3) NUM processors
####  4) out reference genome directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

while getopts a:g:p:d: option
  do
    case "${option}"
      in
      a) GTF=${OPTARG};;
      g) FASTA=${OPTARG};;
      p) THREADS=${OPTARG};;
      d) GENOME_DIR=${OPTARG};;
    esac
done

mkdir $GENOME_DIR/star_index

/mnt/HDD2/colin/bin/STAR-2.7.2d/bin/Linux_x86_64/STAR \
     --runThreadN $THREADS\
     --runMode genomeGenerate \
     --sjdbOverhang 124\
     --genomeChrBinNbits 16 \
     --genomeDir $GENOME_DIR/star_index \
     --genomeFastaFiles $FASTA \
     --sjdbGTFfile $GTF

# END
