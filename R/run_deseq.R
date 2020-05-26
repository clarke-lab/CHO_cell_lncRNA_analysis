#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: run_deseq2.r
##
## Purpose of script: To carry out count based differential expression analyis 
## for the fully StringTie transcriptome assembly. The differential expression results for 
## lncRNAs are parsed and the classififcaiton and orthology information is added 
##
## Author: NIBRT
##
## Date Created: Dec-2020
##
## Email: colin.clarke@nibrt.ie
##
## -----------------------------------------------------------------------------
##
## Info: Genes differentially expressed with +/- 1.5 fold difference, along with a
## BH adjusted p-value and baseMean >=100 are considered significant
##
## -----------------------------------------------------------------------------


print("The error - The query to the BioMart webservice returned an invalid result: 
      biomaRt expected a character string of length 1 - can result from temporary unavailiability of BioMart. 
      Rerun process_rmats.R or use a different host for the useMart function in rmats annotation.R")

# load libraries
suppressMessages(library("DESeq2"))
suppressMessages(library("biomaRt"))
suppressMessages(library("xlsx"))
suppressMessages(library(ggplot2))
suppressMessages(library(refGenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))

# input arguments
args <- commandArgs(TRUE)
count_dir <- args[1] 
results_dir <- args[2]

# volcano plot function
volcano_plot <- function(de_results) {
  de_results <- de_results[de_results$baseMean >= 100, ]
  
  # identify significant
  de_results$threshold_DE <- de_results$padj < 0.05 & 
                             abs(de_results$log2FoldChange) >= 0.5849625
  volcano_data <- cbind(de_results$log2FoldChange, 
                    de_results$padj, 
                    as.logical(de_results$threshold_DE))
  
  # data frame for all required information
  volcano_data <- data.frame(
    ensemblid = rownames(de_results),
    log2FoldChange = de_results$log2FoldChange,
    padj = de_results$padj,
    threshold = de_results$threshold_DE)
  
  volcano_data <- volcano_data %>% dplyr::mutate(pointcolor = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ "#F9B288",
    log2FoldChange < 0 & threshold == 1 ~ "#A2D9F9",
    threshold == 0 ~ "gray"
  ))
  volcano_data <- volcano_data %>% dplyr::mutate(pointclass = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ "Increase at 31C",
    log2FoldChange < 0 & threshold == 1 ~ "Increase at 37C",
    threshold == 0 ~ "Not Significant"
  ))
  volcano_data <- volcano_data %>% dplyr::mutate(pointsize = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ "0.5",
    log2FoldChange < 0 & threshold == 1 ~ "0.5",
    threshold == 0 ~ "0.2"
  ))
  return(volcano_data)
}

modify_syteny_key <- function(syteny_table) {
  ensembl <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "cgcrigri_gene_ensembl",
    host = "uswest.ensembl.org"
  )
  
  for (i in 1:dim(syteny_table)[1]) {
    if (length(grep("MSTRG", syteny_table[i, 1])) > 0) {
      syteny_table[i, 1] <- substr(syteny_table[i, 1], 1, nchar(syteny_table[i, 1]) - 2) #remove the 
      #print("mstrg")
    } else {
      syteny_table[i, 1] <- getBM(
        attributes = "ensembl_gene_id",
        filters = "ensembl_transcript_id", 
        values = syteny_table[i, 1],
        mart = ensembl, uniqueRows = TRUE
      )$ensembl_gene_id
    }
  }
  return(syteny_table)
}

# create object for differential expression
count_file_names <- grep("counts", 
                         list.files(count_dir), 
                         value = T)

cell_type <- c("NTS", "NTS", "NTS", "NTS",
               "TS", "TS", "TS", "TS")

sample_information <- data.frame(
  sampleName = count_file_names,
  fileName = count_file_names,
  condition = cell_type
)

DESeq_data <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_information,
  directory = count_dir,
  design = ~condition
)

colData(DESeq_data)$condition <- factor(colData(DESeq_data)$condition,
                                          levels = c("TS", "NTS")
)

# set 37C as the comparator condition
DESeq_data$condition <- relevel(DESeq_data$condition, "NTS")

# calculate differential expression using the DESeq wrapper function
suppressMessages(DESeq_data <- DESeq(DESeq_data))

# set differential expression criteria
de_results <- results(DESeq_data,
                      lfcThreshold = 0,
                      independentFiltering = T)

# retain significant results
sig_de_results <- subset(
  de_results,
  abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 100
)

sig_de_results <- sig_de_results[
  order(sig_de_results$log2FoldChange, decreasing = T),]

# identify differentially expressed lncRNA genes
# lncrna transcript list
lncrna_transcripts <- read.table("lncrna_annotation/monoexonic_filter/final_lncrna_list.txt", 
                                 header = F, stringsAsFactors = F)$V1

# load gene GTF
suppressMessages(stringtie_txdb <- makeTxDbFromGFF(
  "transcriptome_assembly/non_protein_coding_stringtie.gtf"
))

# collapse lncRNA transcripts to genes 
suppressMessages(lncrna_genes_all <- select(stringtie_txdb,
                                            keys = lncrna_transcripts,
                                            columns = c("TXCHROM", 
                                                        "GENEID", 
                                                        "TXSTRAND", 
                                                        "TXSTART", 
                                                        "TXEND"), 
                                            keytype = "TXNAME"))

# lncrna genes
lncrna_genes <- lncrna_genes_all[!is.na(lncrna_genes_all$GENEID), ]

# extract only lncRNA genes 
all_lncrna_de_results <- de_results[lncrna_genes$GENEID, ]

# create object for volcano
volcano_data <- volcano_plot(all_lncrna_de_results)


volcano_ggplot <- ggplot(volcano_data) +
  geom_point(aes(x = log2FoldChange, 
                 y = -log10(padj), 
                 color = pointcolor, 
                 fill = pointclass), 
             size = 0.1) + 
  scale_colour_manual(
    name = "the colour",
    breaks = as.factor(volcano_data$pointclass),
    values = c("#F9B288", "#A2D9F9", "gray"),
    labels = as.factor(volcano_data$pointclass)
  ) +
  xlab(bquote(~ Log[2] ~ "fold change")) + ylab("-log10 BH adjusted p-value") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed", color = "gray"
  ) + geom_vline(
    xintercept = c(-0.5849625, 0.5849625),
    linetype = "dashed",
    color = "light gray"
  ) +  
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7)
  )

# save
ggsave(paste(results_dir,"/lncrna_volcano.tiff",sep=""), 
       plot = volcano_ggplot, 
       height = 5, width = 5, units = "in", dpi = 600)


# add additional info to DE lncRNA table
# get significant lncrnas 
sig_de_results_lncrna <- sig_de_results[rownames(sig_de_results) %in% lncrna_genes$GENEID, ]

# get mouse orthology information
cho_mouse_synteny <- read.table("lncrna_annotation/liftover/lncrna_cho_to_mouse_ortho.txt", 
                                header = F, 
                                stringsAsFactors = F)

# collapse the results of transcript orthology to genes
cho_mouse_lncrna_geneid <- modify_syteny_key(cho_mouse_synteny)

# annotate the mouse synteny hits
ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = "uswest.ensembl.org")

mouse_lncrnas <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id_version", values = cho_mouse_lncrna_geneid$V2,
  mart = ensembl, uniqueRows = TRUE
)

n_occur <- data.frame(table(mouse_lncrnas$ensembl_gene_id_version))
mouse_lncrnas_duplicated <- mouse_lncrnas[mouse_lncrnas$ensembl_gene_id_version %in% 
                                            as.character(n_occur[n_occur$Freq > 1, ]$Var1), ]
mouse_lncrnas_duplicated[order(mouse_lncrnas_duplicated$external_gene_name), ]
mouse_lncrnas_deduplicated <- unique(mouse_lncrnas[, 1:3])

# add the external gene names to mouse
rownames(mouse_lncrnas_deduplicated) <- mouse_lncrnas_deduplicated$ensembl_gene_id_version
cho_mouse_lncrna_geneid$external_gene_name <- mouse_lncrnas_deduplicated[cho_mouse_lncrna_geneid$V2, 2]
sig_de_results_lncrna$mouse_ensid <- cho_mouse_lncrna_geneid$V2[match(gsub("\\Q|\\E.*","",rownames(sig_de_results_lncrna)),
                                                                      cho_mouse_lncrna_geneid$V1)]
sig_de_results_lncrna$mouse_symbol <- cho_mouse_lncrna_geneid$external_gene_name[match(gsub("\\Q|\\E.*","",rownames(sig_de_results_lncrna)),
                                                                                       cho_mouse_lncrna_geneid$V1)]

# get mouse orthology information
cho_human_synteny <- read.table("lncrna_annotation/liftover/lncrna_cho_to_human_otho.txt", 
                                header = F, 
                                stringsAsFactors = F)

# collapse the results of transcript orthology to genes
cho_human_synteny_mod <- modify_syteny_key(cho_human_synteny)

# annotate the human synteny hits
ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "uswest.ensembl.org")

biomart_human <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id_version", values = cho_human_synteny_mod$V2,
  mart = ensembl, uniqueRows = TRUE
)

n_occur <- data.frame(table(biomart_human$ensembl_gene_id_version))
biomart_human_duplicated <- biomart_human[biomart_human$ensembl_gene_id_version %in% 
                                            as.character(n_occur[n_occur$Freq > 1, ]$Var1), ]
biomart_human_duplicated[order(biomart_human_duplicated$external_gene_name), ]
biomart_human_deduplicated <- unique(biomart_human[, 1:3])

# add the external gene names to human synteny
rownames(biomart_human_deduplicated) <- biomart_human_deduplicated$ensembl_gene_id_version
cho_human_synteny_mod$external_gene_name <- biomart_human_deduplicated[cho_human_synteny_mod$V2, 2]

sig_de_results_lncrna$human_ensid <- cho_human_synteny_mod$V2[match(gsub("\\Q|\\E.*","",rownames(sig_de_results_lncrna)), 
                                                                    cho_human_synteny_mod$V1)]
sig_de_results_lncrna$human_symbol <- cho_human_synteny_mod$external_gene_name[match(gsub("\\Q|\\E.*","",rownames(sig_de_results_lncrna)), 
                                                                                     cho_human_synteny_mod$V1)]


# add the FEELnc classifiation for lncRNA type
lncRNA_classification <- read.csv("lncrna_annotation/classification/final_classification.txt", 
                                  header = T, stringsAsFactors = F, row.names = 1)[, 2:9]

# classification match columns
match_classification <- lncRNA_classification[match(rownames(sig_de_results_lncrna), 
                                                    lncRNA_classification$lncRNA_gene), 8]
match_classification[is.na(match_classification)] <- "intergenic"
match_classification <- lncRNA_classification[match(rownames(sig_de_results_lncrna), 
                                                    lncRNA_classification$lncRNA_gene), 8]
match_classification[is.na(match_classification)] <- "intergenic"
match_partner <- lncRNA_classification[match(rownames(sig_de_results_lncrna), 
                                             lncRNA_classification$lncRNA_gene), 3]
match_partner[is.na(match_partner)] <- "intergenic"
match_distance <- lncRNA_classification[match(rownames(sig_de_results_lncrna), 
                                              lncRNA_classification$lncRNA_gene), 5]
match_distance[is.na(match_distance)] <- "intergenic"


# add the differential expression results for any protein coding genes published 
# tzani et al 2020
de_protein_coding <- read.xlsx("data/bit27365-sup-0007-Table_S3.xlsx", sheetIndex = 1,header = T)
de_protein_coding <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), ]
de_protein_coding_entrezid <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), 2]
de_protein_coding_symbol <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), 3]
de_protein_coding_description <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), 4]
de_protein_coding_fc <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), 6]
de_protein_coding_padj <- de_protein_coding[match(match_partner, de_protein_coding$Ensembl.Gene.ID), 10]


# create a summary table

lncrna_de_matrix <- data.frame(cbind(
  sig_de_results_lncrna[, c(1:2, 6, 7, 8, 9, 10)],
  match_classification, 
  match_distance,
  match_partner,
  de_protein_coding_symbol,
  de_protein_coding_entrezid,
  de_protein_coding_description,
  de_protein_coding_fc,
  de_protein_coding_padj
))

colnames(lncrna_de_matrix) <- c(
  "lncRNA baseMean",
  "lncRNA log2FC",
  "lncRNA BH Adjusted P-value",
  "Orthologus Mm10 ENSGID",
  "Orthologus Mm10 Symbol",
  "Orthologus Hg38 ENSGID",
  "Orthologus Hg38 Symbol",
  "lncRNA type",
  "Match distance",
  "CHOK1 ENSGID",
  "CHOK1 Symbol",
  "CHOK1 ENTREZID",
  "CHOK1 Description",
  "mRNA log2FC",
  "lncRNA BD Adjusted P-value"
)


# save the de_results object
#save(de_results, file = paste(results_dir, "/de_results.rData", sep=""))
fn <- paste(results_dir, "/Annotated_lncRNA_DESeq2_results.xlsx",sep="")
suppressMessages(if (file.exists(fn)) {file.remove(fn)})

# save the results to Excel
write.xlsx(lncrna_de_matrix,
           file = fn,
           sheetName = "DE lncRNAs", append = TRUE
)
