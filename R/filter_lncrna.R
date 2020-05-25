#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: filter_lncrna.R
##
## Purpose of script: compare the raw feelnc results with other metrics of coding probability
##
## Author: NIBRT
##
## Date Created: Dec-2020
##
## Email: colin.clarke@nibrt.ie
##
## -----------------------------------------------------------------------------
##
## Info: the following code filters the feelnc identified transcripts following comparison with
## PSI >= 0.1 and an FDR < 0.05 to be considered signifciant
## Each rMats results is annotated using bioMart, to resolve transcripts with sense
## overlap the exon coordinates are searched
## ---------------------------------------------------------------------------
suppressMessages(library(refGenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))

# command line arguments
args <- commandArgs(TRUE)
candidate_lncrna_gtf <- args[1]
cpc_output <- args[2]
cpat_output <- args[3]
swissprot_blast_output <- args[4]
mirbase_blast_output <- args[5]
pfam_domain_hit_output <- args[6]
rfam_blast_output <- args[7]
tpm_expression_matrix <- args[8]
feelnc_lncrna_classification <- args[9]
stringtie_gtf <- args[10]
output_dir <- args[11]

# make a folder to hold the results of this filter
system(paste("mkdir ", output_dir, "/firstpass_filter", sep = ""))

# load the feelnc gtf
feelnc_gtf <- rtracklayer::import(candidate_lncrna_gtf)

# remove mtDNA transcripts
feelnc_gtf <- feelnc_gtf[seqnames(feelnc_gtf) != "MT"]
feelnc_lncrna_transcripts <- unique(feelnc_gtf$transcript_id)

write("---------------FEELNC-------------------", stdout())
write(
  paste(length(feelnc_lncrna_transcripts),
    "candidate lncRNAs identified by FEELNc",
    sep = " "
  ),
  stdout()
)

# import the lncRNAs annotated in Ensembl 
ensembl_lncrna_transcripts <- read.table("reference_genome/ensembl_lncrna_transcript.txt", stringsAsFactors = F)$V1
in_ensembl <- intersect(feelnc_lncrna_transcripts, ensembl_lncrna_transcripts)

write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

# coding potential calculator
# CPC2 results
cpc2 <- read.table(cpc_output, sep = "\t")
cpc2 <- cpc2[cpc2$V8 == "noncoding", ]
cpc2_lncrna_transcripts <- as.character(cpc2$V1)

# CPAT  define coding/non-coding based on mouse (0.44)
cpat <- read.table(cpat_output)
cpat <- cpat[cpat$coding_prob < 0.44, ]
cpat_lncrna_transcripts <- rownames(cpat)

# overlap results of 3 coding potential calculators
lncrna_all_noncoding <- intersect(
  intersect(
    feelnc_lncrna_transcripts,
    cpc2_lncrna_transcripts
  ),
  cpat_lncrna_transcripts
)

write("---------------CPAT,CPC2-------------------", stdout())
write(
  paste(length(lncrna_all_noncoding),
    "candidate lncRNAs common to FEELNc, CPC2 & CPAT",
    sep = " "
  ),
  stdout()
)

in_ensembl <- intersect(lncrna_all_noncoding, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

# filter results with a blast hit against SWISSPROT
blastp <- read.table(swissprot_blast_output, sep = "\t")
protein_hit <- data.frame("data" = as.character(blastp$V1))
protein_hit <- substr(
  as.character(protein_hit$data),
  1, nchar(as.character(protein_hit$data)) - 3
)
no_blastp_hit <- !(lncrna_all_noncoding %in% protein_hit)
lncrna_all_noncoding <- lncrna_all_noncoding[no_blastp_hit]

write("---------------SWISSPROT BLAST-------------------", stdout())
write(
  paste(length(lncrna_all_noncoding),
    "candidate lncRNAs remaining after SWISS prot filtering",
    sep = " "
  ),
  stdout()
)

in_ensembl <- intersect(lncrna_all_noncoding, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)


# filter results with blast hit against miRBase
blastn_mir <- read.table(mirbase_blast_output, sep = "\t")
mirna_hit <- as.character(blastn_mir$V1)

no_mirbase_hit <- !(lncrna_all_noncoding %in% mirna_hit)
lncrna_all_noncoding <- lncrna_all_noncoding[no_mirbase_hit]
write("---------------MIRBASE BLAST-------------------", stdout())
write(
  paste(length(lncrna_all_noncoding),
    "candidate lncRNAs remaining after mirbase filtering",
    sep = " "
  ),
  stdout()
)

in_ensembl <- intersect(lncrna_all_noncoding, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

# filter results with PFAM domain
pfam_domain_hit <- read.table(pfam_domain_hit_output)
pfam_domain_hit <- as.character(pfam_domain_hit[, 1])
no_pfam_domain_hit <- !(lncrna_all_noncoding %in% pfam_domain_hit)
lncrna_all_noncoding <- lncrna_all_noncoding[no_pfam_domain_hit]

write("---------------PFAM-------------------", stdout())
write(
  paste(length(lncrna_all_noncoding),
    "candidate lncRNAs remaining after PFAM domain filtering",
    sep = " "
  ),
  stdout()
)

in_ensembl <- intersect(lncrna_all_noncoding, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

# filter results with RFAM(-lncRNAs sequences)
blastn_rfam <- read.table(rfam_blast_output, sep = "\t")
ncrna_hit <- as.character(blastn_rfam$V1)
no_rfam_hit <- !(lncrna_all_noncoding %in% ncrna_hit)
# intersect(lncrna_all_noncoding,ncrna_hit)
lncrna_all_noncoding <- lncrna_all_noncoding[no_rfam_hit]

write("---------------RFAM-------------------", stdout())
write(
  paste(length(lncrna_all_noncoding),
    "candidate lncRNAs remaining after RFAM filtering",
    sep = " "
  ),
  stdout()
)

in_ensembl <- intersect(lncrna_all_noncoding, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2974 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

# expression filter
tpm_matrix <- read.table(tpm_expression_matrix, sep = "\t", header = T, row.names = 1)
lncrna_tpm_matrix <- tpm_matrix[lncrna_all_noncoding, ]
tpm_filter <- apply(lncrna_tpm_matrix, 1, max) > 1 # greater than TPM in 1 sample
lncrna_transcripts <- rownames(lncrna_tpm_matrix[tpm_filter, ])
write("---------------TPM FILTER-------------------", stdout())
write(
  paste(length(lncrna_transcripts),
    "candidate lncRNAs with TPM > 1 in at least 1 sample",
    sep = " "
  ),
  stdout()
)
in_ensembl <- intersect(lncrna_transcripts, ensembl_lncrna_transcripts)
write(
  paste(length(in_ensembl),
    "of 2975 lncRNAs annotated by ENSEMBL present ",
    sep = " "
  ),
  stdout()
)

write("-----------------------------------------", stdout())

suppressMessages(stringtie_txdb <- makeTxDbFromGFF(
  stringtie_gtf
))

suppressMessages(lncrna_genes <- select(stringtie_txdb,
  keys = lncrna_transcripts,
  columns = "GENEID", keytype = "TXNAME"
))

save(lncrna_transcripts, lncrna_genes, file = paste(output_dir, "firstpass_filter/lncRNA_filtering.rData", sep = "/"))

write(setdiff(subset(lncrna_transcripts, grepl("ENSCGRT", lncrna_transcripts)), in_ensembl), paste(output_dir, "/firstpass_filter/ensembl_transcripts_other_rna.txt", sep = ""))

write("writing novel and complete lists of files", stdout())
write(subset(lncrna_transcripts, !grepl("ENSCGRT", lncrna_transcripts)), paste(output_dir, "/firstpass_filter/novel_lncrna_transcripts.txt", sep = ""))
write(lncrna_transcripts, paste(output_dir, "/firstpass_filter/all_lncrna_transcripts.txt", sep = ""))


system(paste("grep -wFf ", output_dir, "/firstpass_filter/all_lncrna_transcripts.txt", " transcriptome_assembly/non_protein_coding_stringtie.gtf >", output_dir, "/firstpass_filter/all_lncrna.gtf", sep = ""))
# system(paste("grep -wFf reference_genome/ensembl_lncrna_transcript.txt", " stringtie_output/non_protein_coding_stringtie.gtf >", output_dir, "/firstpass_filter/ensembl_lncrna_non_protein_coding_stringtie.gtf" ,sep=""))

# use transdecoder to create a fasta file of all lncRNAs
system(paste("~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl ", output_dir, "/firstpass_filter/all_lncrna.gtf reference_genome/ensembl_chok1_genome.fa > ", output_dir, "/firstpass_filter/lncRNA.fasta", sep = ""))

write("-----------------------------------------", stdout())
write("complete: Data saved", stdout())
write(paste("All transcripts:", output_dir, "/firstpass_filter/all_lncrna_transcripts.txt", sep = ""), stdout())
write(paste("Novel transcripts:", output_dir, "/firstpass_filter/novel_lncrna_transcripts.txt", sep = ""), stdout())
write(paste("R data object:", output_dir, "/firstpass_filter/lncRNA_filtering.rData", sep = ""), stdout())
write(paste("Ensembl RNA transcripts annotated as lncRNA:", output_dir, "/firstpass_filter/ensembl_transcripts_other_rna.txt", sep = ""), stdout())
write(paste("lncRNA transcripts sequences:", output_dir, "/firstpass_filter/lncRNA.fasta", sep = ""), stdout())
write("-------------Summary-------------------------", stdout())
write(
  paste(length(lncrna_transcripts),
    " lncRNA transcripts were identified",
    sep = " "
  ),
  stdout()
)
write(
  paste(length(subset(lncrna_transcripts, !grepl("ENSCGRT", lncrna_transcripts))),
    " novel lncRNA transcripts were identified",
    sep = " "
  ),
  stdout()
)

write(
  paste(length(in_ensembl),
    "of lncRNAs previously annotated by ENSEMBL",
    sep = " "
  ),
  stdout()
)
write(paste(length(setdiff(subset(lncrna_transcripts, grepl("ENSCGRT", lncrna_transcripts)), in_ensembl)),
  " other RNAs classified as lncRNA",
  sep = ""
), stdout())

write("-----------------------------------------", stdout())
write(paste(length(feelnc_lncrna_transcripts) - length(intersect(feelnc_lncrna_transcripts, cpc2_lncrna_transcripts))), "filtered by cpc2")
write(paste(length(feelnc_lncrna_transcripts) - length(intersect(feelnc_lncrna_transcripts, cpat_lncrna_transcripts))), "filtered by cpat")
# write(paste(length(feelnc_lncrna_transcripts)-length(intersect(feelnc_lncrna_transcripts,cpat_lncrna_transcripts))), "filtered by cpat")
# write(paste(length(feelnc_lncrna_transcripts)-length(intersect(feelnc_lncrna_transcripts,cpat_lncrna_transcripts))), "filtered by cpat")
