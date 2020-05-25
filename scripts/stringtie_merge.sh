#!/usr/bin/env bash
#### Merge the individual transcript assemblies
#### inputs are: 1) assembled transcript parent directory and 2) the refernece GTF
#### 3) reference genome directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-t = assembled transcript directory"
        echo "-g = path to reference annotation"
        echo "-r = reference genome directory"
        exit 2
fi
while getopts t:g:r: option
  do
    case "${option}"
      in
      t) TRANSCRIPT_DIR=${OPTARG};;
      g) GTF=${OPTARG};;
      r) REF_DIR=${OPTARG};;
    esac
done


readlink -f $TRANSCRIPT_DIR/individual_gtfs/*.gtf > $TRANSCRIPT_DIR/mergelist.txt

/mnt/HDD2/colin/bin/stringtie-2.0.3.Linux_x86_64/stringtie \
--merge $TRANSCRIPT_DIR/mergelist.txt \
-o $TRANSCRIPT_DIR/stringtie_original.gtf \
-G $GTF \
-f 0.1  \
-c 10

# create a file liniking stringtie ID to ENSEMBL geNE ID
grep -wFf $REF_DIR/protein.coding.genes.txt $TRANSCRIPT_DIR/stringtie_original.gtf | \
grep -v exon | awk '{print $10, $NF}' | uniq | tr -d \" | tr -d \; > \
$TRANSCRIPT_DIR/stringtie_ensembl_gene_mapping.txt


# make a version of the stringtie assembly for differential exression containing
# transcripts
perl scripts/mstrg_prep.pl $TRANSCRIPT_DIR/stringtie_original.gtf > $TRANSCRIPT_DIR/stringtie_original_appended.gtf

grep 'MSTRG.*|ENSCGRG.*|ENSC.*' $TRANSCRIPT_DIR/stringtie_original_appended.gtf | \
grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt

grep -v -F -f $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt $TRANSCRIPT_DIR/stringtie_original_appended.gtf > \
$TRANSCRIPT_DIR/stringtie_original.appended.fp.filtered.gtf

# remove transcripts without strand information
awk '$7 != "." {print}' $TRANSCRIPT_DIR/stringtie_original.appended.fp.filtered.gtf > \
$TRANSCRIPT_DIR/stringtie.all.transcripts.gtf

# create a stringtie GTF of protein coding genes only

grep -wFf $REF_DIR/protein.coding.genes.txt $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf > \
$TRANSCRIPT_DIR/stringtie.protein.coding.gtf

# remove intermediate files
rm $TRANSCRIPT_DIR/stringtie_original_appended.gtf \
$TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt \
$TRANSCRIPT_DIR/stringtie_original.appended.fp.filtered.gtf \

# make to GTF file to search for lncRNAs
# filter transcripts with 1bp exon overlap with annotated protein coding genes
#create bedtools overlap, only for lncRNA
awk '{if($3=="exon"){print $10}}' $TRANSCRIPT_DIR/stringtie_original.gtf | \
sed 's/"//g;s/;//g' | sort | uniq > $TRANSCRIPT_DIR/stringtie.gene.txt

awk '{if($3=="exon"){print}}' $TRANSCRIPT_DIR/stringtie_original.gtf > \
$TRANSCRIPT_DIR/stringtie.exon.gtf

grep "protein_coding" $GTF | awk '{if($3=="exon"){print}}' > \
$TRANSCRIPT_DIR/coding.exon.gtf

bedtools intersect -s -u -a $TRANSCRIPT_DIR/stringtie.exon.gtf \
                         -b $TRANSCRIPT_DIR/coding.exon.gtf > \
                          $TRANSCRIPT_DIR/overlap_exon.gtf

awk '{print $10}' $TRANSCRIPT_DIR/overlap_exon.gtf | sed 's/"//g;s/;//g' | sort | uniq  > \
$TRANSCRIPT_DIR/overlap.exon.gene.id.txt

comm -23 $TRANSCRIPT_DIR/stringtie.gene.txt $TRANSCRIPT_DIR/overlap.exon.gene.id.txt | \
grep -wFf - $TRANSCRIPT_DIR/stringtie_original.gtf > $TRANSCRIPT_DIR/strigtie_overlap_filtered.gtf

# append ensembl gene ids to MSTRG GTF
perl scripts/mstrg_prep.pl $TRANSCRIPT_DIR/strigtie_overlap_filtered.gtf > \
$TRANSCRIPT_DIR/stringtie_merged.appended.gtf

# find instances where stringtie has asemebled transcripts from 2 or more overlaping loci and created a new "gene".
# The final field of the GTF file will contain an MSTRG ID not an ENS ID
grep 'MSTRG.*|ENSCGRG.*|ENSC.*' $TRANSCRIPT_DIR/stringtie_merged.appended.gtf | \
grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt

# remove assembled transcripts spanning two or more sense overlapping genes transcripts
grep -v -F -f $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt $TRANSCRIPT_DIR/stringtie_merged.appended.gtf > \
$TRANSCRIPT_DIR/stringtie_merged.appended.fp.filtered.gtf

# remove transcripts without strand
awk '$7 != "." {print}' $TRANSCRIPT_DIR/stringtie_merged.appended.fp.filtered.gtf > $TRANSCRIPT_DIR/non_protein_coding_stringtie.gtf

# compare to the CHO K1 reference sequence
gffcompare \
-o $TRANSCRIPT_DIR/gffcompare \
-r $GTF $TRANSCRIPT_DIR/non_protein_coding_stringtie.gtf

# END
