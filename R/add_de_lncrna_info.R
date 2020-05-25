## subset and annotate DESeq results
suppressMessages(library(biomaRt))
suppressMessages(library(refGenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library("DESeq2"))

modify_syteny_key <- function(syteny_table) {
  ensembl <- useMart("ensembl")

  ensembl <- useDataset("cgcrigri_gene_ensembl", mart = ensembl)
  filters <- listFilters(ensembl)
  attributes <- listAttributes(ensembl)

  for (i in 1:dim(syteny_table)[1]) {
    if (length(grep("MSTRG", syteny_table[i, 1])) > 0) {
      syteny_table[i, 1] <- substr(syteny_table[i, 1], 1, nchar(syteny_table[i, 1]) - 2)
      print("mstrg")
    } else {
      syteny_table[i, 1] <- getBM(
        attributes = "ensembl_gene_id",
        filters = "ensembl_transcript_id", values = syteny_table[i, 1],
        mart = ensembl, uniqueRows = TRUE
      )$ensembl_gene_id
      print("complete")
    }
  }
  return(syteny_table)
}


volcano_plot <- function(de_results) {
  de_results <- de_results[de_results$baseMean >= 100, ]

  # add significant
  de_results$threshold_DE <- de_results$padj < 0.05 & abs(de_results$log2FoldChange) >= 0.5849625

  DE_event <- cbind(de_results$log2FoldChange, de_results$padj, as.logical(de_results$threshold_DE))

  # data frame for all required information
  DE_event <- data.frame(
    ensemblid = rownames(de_results),
    log2FoldChange = de_results$log2FoldChange,
    padj = de_results$padj,
    threshold = de_results$threshold_DE
  )

  DE_event <- DE_event %>% mutate(pointcolor = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "#F9B288",
    log2FoldChange < 0 & threshold == 1 ~ "#A2D9F9",
    threshold == 0 ~ "gray"
  ))
  DE_event <- DE_event %>% mutate(pointclass = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "Increase at 31C",
    log2FoldChange < 0 & threshold == 1 ~ "Increase at 37C",
    threshold == 0 ~ "Not Significant"
  ))
  DE_event <- DE_event %>% mutate(pointsize = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "0.5",
    log2FoldChange < 0 & threshold == 1 ~ "0.5",
    threshold == 0 ~ "0.2"
  ))

  return(DE_event)
}

lncrna_genes$GENEID <- lncrna_genes$GENEID[is.na(lncrna_genes$GENEID), ]


# load DESeq results
load("differential_expression/DESeq2_results/de_results.RData")
# lncRNA annotation info
lncrna_transcripts <- read.table("lncrna_annotation/monoexonic_filter/final_lncrna_list.txt", header = F, stringsAsFactors = F)$V1

suppressMessages(stringtie_txdb <- makeTxDbFromGFF(
  "stringtie/non_protein_coding_stringtie.gtf"
))

suppressMessages(lncrna_genes_all <- select(stringtie_txdb,
  keys = lncrna_transcripts,
  columns = c("TXCHROM", "GENEID", "TXSTRAND", "TXSTART", "TXEND"), keytype = "TXNAME"
))

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
lncrna_genes <- lncrna_genes_all[!is.na(lncrna_genes_all$GENEID), ]

# volcano plotting
all_lncrna_de_results <- de_results[lncrna_genes$GENEID, ]

DE_event <- volcano_plot(all_lncrna_de_results)


volcano_ggplot <- ggplot(DE_event) +
  geom_point(aes(x = log2FoldChange, 
                 y = -log10(padj), 
                 color = pointcolor, 
                 fill = pointclass), 
             size = 0.1) + 
  scale_colour_manual(
    name = "the colour",
    breaks = as.factor(DE_event$pointclass),
    values = c("#F9B288", "#A2D9F9", "gray"),
    labels = as.factor(DE_event$pointclass)
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
  
  
  

volcano_ggplot <- volcano_ggplot + scale_colour_manual(
  name = "the colour",
  breaks = as.factor(DE_event$pointclass),
  values = c("#F9B288", "#A2D9F9", "gray"),
  labels = as.factor(DE_event$pointclass)
)

volcano_ggplot <- volcano_ggplot + xlab(bquote(~ Log[2] ~ "fold change")) + ylab("-log10 BH adjusted p-value")
volcano_ggplot <- volcano_ggplot + geom_hline(
  yintercept = -log10(0.05),
  linetype = "dashed", color = "gray"
)
volcano_ggplot <- volcano_ggplot + geom_vline(
  xintercept = c(-0.5849625, 0.5849625),
  linetype = "dashed",
  color = "light gray"
)
volcano_ggplot <- volcano_ggplot + theme(
  legend.position = "none",
  plot.title = element_text(size = rel(1.5), hjust = 0.5),
  axis.title = element_text(size = rel(1.25))
)

volcano_ggplot <- volcano_ggplot + theme_bw()
volcano_ggplot <- volcano_ggplot + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 7),
  legend.position = "none"
)

ggsave(paste("differential_expression/volcano.tiff"), plot = volcano_ggplot, height = 5, width = 5, units = "in", dpi = 600)



# select annotated lncrnas for significant de transcripts
sig_de_results_lncrna <- sig_de_results[rownames(sig_de_results) %in% lncrna_genes$GENEID, ]

# required files for annotation
cho_mouse_synteny <- read.table("lncrna_annotation/liftover/lncrna_cho_to_mouse_ortho.txt", header = F, stringsAsFactors = F)
cho_human_synteny <- read.table("lncrna_annotation/liftover/lncrna_cho_to_human_otho.txt", header = F, stringsAsFactors = F)

cho_mouse_synteny_mod <- modify_syteny_key(cho_mouse_synteny)
cho_human_synteny_mod <- modify_syteny_key(cho_human_synteny)

# annotate the mouse synteny hits
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

biomart_mouse <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id_version", values = cho_mouse_synteny_mod$V2,
  mart = ensembl, uniqueRows = TRUE
)

n_occur <- data.frame(table(biomart_mouse$ensembl_gene_id_version))
biomart_mouse_duplicated <- biomart_mouse[biomart_mouse$ensembl_gene_id_version %in% as.character(n_occur[n_occur$Freq > 1, ]$Var1), ]
biomart_mouse_duplicated[order(biomart_mouse_duplicated$external_gene_name), ]
biomart_mouse_deduplicated <- unique(biomart_mouse[, 1:3])

# add the external gene names to mouse
rownames(biomart_mouse_deduplicated) <- biomart_mouse_deduplicated$ensembl_gene_id_version
cho_mouse_synteny_mod$external_gene_name <- biomart_mouse_deduplicated[cho_mouse_synteny_mod$V2, 2]
sig_de_results_lncrna$mouse_ensid <- cho_mouse_synteny_mod$V2[match(rownames(sig_de_results_lncrna), cho_mouse_synteny_mod$V1)]
sig_de_results_lncrna$mouse_symbol <- cho_mouse_synteny_mod$external_gene_name[match(rownames(sig_de_results_lncrna), cho_mouse_synteny_mod$V1)]

# annotate the human synteny hits
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

biomart_human <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id_version", values = cho_human_synteny_mod$V2,
  mart = ensembl, uniqueRows = TRUE
)

n_occur <- data.frame(table(biomart_human$ensembl_gene_id_version))
biomart_human_duplicated <- biomart_human[biomart_human$ensembl_gene_id_version %in% as.character(n_occur[n_occur$Freq > 1, ]$Var1), ]
biomart_human_duplicated[order(biomart_human_duplicated$external_gene_name), ]
biomart_human_deduplicated <- unique(biomart_human[, 1:3])

# add the external gene names to mouse synteny
rownames(biomart_human_deduplicated) <- biomart_human_deduplicated$ensembl_gene_id_version
cho_human_synteny_mod$external_gene_name <- biomart_human_deduplicated[cho_human_synteny_mod$V2, 2]

sig_de_results_lncrna$human_ensid <- cho_human_synteny_mod$V2[match(rownames(sig_de_results_lncrna), cho_human_synteny_mod$V1)]
sig_de_results_lncrna$human_symbol <- cho_human_synteny_mod$external_gene_name[match(rownames(sig_de_results_lncrna), cho_human_synteny_mod$V1)]

lncRNA_classification <- read.csv("lncrna_annotation/classification/final_classification.txt", header = T, stringsAsFactors = F, row.names = 1)[, 2:9]

# classification match columns
match_classification <- lncRNA_classification[match(rownames(sig_de_results_lncrna), lncRNA_classification$lncRNA_gene), 8]
match_classification[is.na(match_classification)] <- "intergenic"

match_classification <- lncRNA_classification[match(rownames(sig_de_results_lncrna), lncRNA_classification$lncRNA_gene), 8]
match_classification[is.na(match_classification)] <- "intergenic"
match_partner <- lncRNA_classification[match(rownames(sig_de_results_lncrna), lncRNA_classification$lncRNA_gene), 3]
match_partner[is.na(match_partner)] <- "intergenic"
match_distance <- lncRNA_classification[match(rownames(sig_de_results_lncrna), lncRNA_classification$lncRNA_gene), 5]
match_distance[is.na(match_distance)] <- "intergenic"

de_protein_coding <- read.csv("differential_expression/31_v_37.csv", header = T, stringsAsFactors = F, row.names = 1)
de_protein_coding <- de_protein_coding[match(match_partner, DESeq2_analysis$...1rownames(de_protein_coding)), ]

de_protein_coding_entrezid <- de_protein_coding[match(match_partner, rownames(de_protein_coding)), 1]
de_protein_coding_symbol <- de_protein_coding[match(match_partner, rownames(de_protein_coding)), 2]
de_protein_coding_description <- de_protein_coding[match(match_partner, rownames(de_protein_coding)), 3]
de_protein_coding_fc <- de_protein_coding[match(match_partner, rownames(de_protein_coding)), 6]
de_protein_coding_padj <- de_protein_coding[match(match_partner, rownames(de_protein_coding)), 10]

# create a summary table

lncrna_de_matrix <- data.frame(cbind(
  sig_de_results_lncrna[, c(1:2, 6, 7, 10, 8, 9)],
  match_classification, match_distance,
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
  "lncRNA BD Adjusted P-value",
  "Orthologus Mm10 ENSGID",
  "Orthologus Mm10 Symbol",
  "Orthologus Hg38 ENSGID",
  "Orthologus Hg38 Symbol",
  "lncRNA type",
  "CHOK1 ENSGID",
  "CHOK1 Symbol",
  "CHOK1 ENTREZID",
  "CHOK1 Description",
  "mRNA log2FC",
  "lncRNA BD Adjusted P-value"
)

write.csv(lncrna_de_matrix, file = "differential_expression/lncrna_differential_expression.csv")
