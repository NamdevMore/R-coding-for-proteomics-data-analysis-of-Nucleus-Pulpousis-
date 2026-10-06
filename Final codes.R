# Set working directory to your project folder
path <- "/Users/moren/Library/CloudStorage/OneDrive-Cedars-SinaiHealthSystem/Sheyn, Dima's files - Box - Sheyn Lab/Projects/Ongoing projects/Microgels - FibMA/Acute stress proteomics code"
setwd(path)

# Confirm it worked
getwd()




pkgs <- c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "matrixStats")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(matrixStats)

# Import the file
raw <- read.csv("report.pg_matrix_GK.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Check dimensions and column names
dim(raw)
colnames(raw)


gene_col <- "Genes"

sample_cols <- c("NPC1","NPC2","NPC3",
                 "sNPC1","sNPC2","sNPC3",
                 "NPC-MG1","NPC-MG2","NPC-MG3",
                 "sNPC-MG1","sNPC-MG2","sNPC-MG3")

# Confirm all sample columns exist in raw
all(sample_cols %in% colnames(raw))

# Build metadata table
meta <- data.frame(Sample = sample_cols) %>%
  mutate(Group = case_when(
    grepl("^sNPC-MG", Sample) ~ "sNPC-MG",
    grepl("^NPC-MG",  Sample) ~ "NPC-MG",
    grepl("^sNPC",    Sample) ~ "sNPC",
    grepl("^NPC",     Sample) ~ "NPC",
    TRUE ~ NA_character_
  ))
meta$Group <- factor(meta$Group, levels = c("NPC", "sNPC", "NPC-MG", "sNPC-MG"))
rownames(meta) <- meta$Sample

print(meta)


# Aggregate duplicate gene entries (average if a gene appears more than once)
expr_agg <- raw %>%
  filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
  select(all_of(c(gene_col, sample_cols))) %>%
  group_by(across(all_of(gene_col))) %>%
  summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

# Build numeric matrix with gene names as rownames
expr_mat <- as.matrix(expr_agg[, -1])
rownames(expr_mat) <- expr_agg[[gene_col]]

# Check dimensions and missingness
dim(expr_mat)
mean(is.na(expr_mat) | is.nan(expr_mat))


# Log2 transform
log_mat <- log2(expr_mat)
log_mat[is.infinite(log_mat)] <- NA

# Check missingness per protein (row)
missing_frac <- rowMeans(is.na(log_mat))
summary(missing_frac)

# Filter: keep proteins detected in at least ~70% of samples (i.e., <=30% missing)
log_mat_filt <- log_mat[missing_frac <= 0.3, ]
dim(log_mat_filt)


# Impute remaining NAs with a low value (simulates "below detection limit")
min_val <- min(log_mat_filt, na.rm = TRUE)

log_mat_imp <- log_mat_filt
log_mat_imp[is.na(log_mat_imp)] <- min_val - 1

# Remove any zero-variance rows (would break heatmap scaling)
row_sd <- apply(log_mat_imp, 1, sd)
log_mat_imp <- log_mat_imp[row_sd > 0, ]

dim(log_mat_imp)


library(pheatmap)
library(RColorBrewer)

# Annotation for columns (samples) showing group membership
annotation_col <- data.frame(Group = meta$Group)
rownames(annotation_col) <- meta$Sample

# Consistent colors for each group
group_colors <- c("NPC" = "#F8766D", "sNPC" = "#7CAE00",
                  "NPC-MG" = "#00BFC4", "sNPC-MG" = "#C77CFF")
ann_colors <- list(Group = group_colors)

# Select top 50 most variable proteins for a readable heatmap
n_top <- 50
gene_var <- apply(log_mat_imp, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:n_top]
top_mat <- log_mat_imp[top_genes, ]

# Plot heatmap: each column = one sample, top annotation = group
pheatmap(top_mat,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,      # keep samples grouped in original order
         cluster_rows = TRUE,       # cluster proteins by expression pattern
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 6,
         main = paste("Top", n_top, "Most Variable Proteins (per sample)"))



# Compute group means for each protein
group_means <- sapply(levels(meta$Group), function(g) {
  cols <- meta$Sample[meta$Group == g]
  rowMeans(log_mat_imp[, cols, drop = FALSE])
})

# group_means is now a matrix: rows = proteins, columns = 4 groups
dim(group_means)
head(group_means)



# Select top most variable proteins across group means
n_top_grp <- 50
grp_var <- apply(group_means, 1, var)
top_genes_grp <- names(sort(grp_var, decreasing = TRUE))[1:n_top_grp]
top_mat_grp <- group_means[top_genes_grp, ]

# Annotation for columns (groups) - reuse same colors as before
annotation_col_grp <- data.frame(Group = factor(colnames(top_mat_grp),
                                                levels = c("NPC","sNPC","NPC-MG","sNPC-MG")))
rownames(annotation_col_grp) <- colnames(top_mat_grp)

pheatmap(top_mat_grp,
         scale = "row",
         annotation_col = annotation_col_grp,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,     # keep NPC -> sNPC -> NPC-MG -> sNPC-MG order
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 6,
         main = paste("Top", n_top_grp, "Variable Proteins (Group Averages)"))


pheatmap(top_mat_grp,
         scale = "row",
         annotation_col = annotation_col_grp,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize = 12,          # overall base font size
         fontsize_row = 10,      # gene name labels
         fontsize_col = 12,      # group name labels
         main = paste("Top", n_top_grp, "Variable Proteins (Group Averages)"))



library(grid)

p <- pheatmap(top_mat_grp,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 12,
              main = paste("Top", n_top_grp, "Variable Proteins (Group Averages)"),
              silent = TRUE)

# Make row and column labels bold
p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 12)

grid.newpage()
grid.draw(p$gtable)


p <- pheatmap(top_mat_grp,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 12,
              cellheight = 12,     # was 9 -> small bump for a bit more row spacing
              cellwidth = 18,
              border_color = "grey60",
              main = paste("Top", n_top_grp, "Variable Proteins (Group Averages)"),
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 7)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)



# Select top 25 most variable proteins across group means
n_top_grp <- 25
grp_var <- apply(group_means, 1, var)
top_genes_grp <- names(sort(grp_var, decreasing = TRUE))[1:n_top_grp]
top_mat_grp <- group_means[top_genes_grp, ]

# Plot heatmap with same publication-style formatting
p <- pheatmap(top_mat_grp,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 11,
              cellheight = 14,       # slightly bigger cells since fewer rows
              cellwidth = 70,
              treeheight_row = 30,
              border_color = "white",
              main = paste("Top", n_top_grp, "Variable Proteins (Group Averages)"),
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

# Recalculate PDF dimensions to match
n_rows <- nrow(top_mat_grp)
n_cols <- ncol(top_mat_grp)

pdf_height <- (n_rows * 14) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("group_heatmap_top25.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()


# Define a curated set of core ECM / matrisome-associated genes
# (expand this list based on your species/tissue if needed)
ecm_genes <- c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL5A1","COL6A1",
               "FN1","LAMA1","LAMB1","LAMC1","ELN","FBLN5","FBLN1","FBN1",
               "SPARC","POSTN","THBS1","THBS2","VTN","TNC","BGN","DCN",
               "LUM","VCAN","HSPG2","NID1","NID2","MMP2","MMP9","TIMP1",
               "TIMP2","PCOLCE","PCOLCE2","CILP","CILP2","SERPINE2","TGM2",
               "TGM3","FGA","FGB","FGG","VWF","CLU","F13B","PLVAP",
               "SELENOP","KRT71","MEAP4")

# Intersect with genes actually present in your data
ecm_present <- intersect(ecm_genes, rownames(group_means))

cat("ECM genes found in dataset:", length(ecm_present), "\n")
print(ecm_present)

# Rank these by variance across groups and take top 10
ecm_var <- apply(group_means[ecm_present, , drop = FALSE], 1, var)
top10_ecm <- names(sort(ecm_var, decreasing = TRUE))[1:min(10, length(ecm_present))]

top_mat_ecm <- group_means[top10_ecm, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_ecm,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 11,
              cellheight = 20,
              cellwidth = 70,
              treeheight_row = 25,
              border_color = "white",
              main = "Top 10 ECM Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_ecm)
n_cols <- ncol(top_mat_ecm)
pdf_height <- (n_rows * 20) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("ecm_top10_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()


# Define a curated set of core ECM / matrisome-associated genes
# (expand/replace with your own list or the Matrisome/Naba reference set for accuracy)
ecm_genes <- c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL5A1","COL6A1",
               "COL6A2","COL6A3","COL7A1","COL12A1","COL14A1",
               "FN1","LAMA1","LAMB1","LAMC1","LAMA2","LAMA4","LAMB2",
               "ELN","FBLN5","FBLN1","FBN1","FBN2",
               "SPARC","POSTN","THBS1","THBS2","THBS4","VTN","TNC",
               "BGN","DCN","LUM","VCAN","HSPG2","NID1","NID2",
               "MMP1","MMP2","MMP3","MMP9","MMP14","TIMP1","TIMP2","TIMP3",
               "PCOLCE","PCOLCE2","CILP","CILP2","SERPINE2","SERPINH1",
               "TGM2","TGM3","FGA","FGB","FGG","VWF","CLU","F13B",
               "PLVAP","SELENOP","KRT71","MEAP4","LOXL1","LOXL2","P4HA1","P4HB")

# Intersect with genes actually present in your data
ecm_present <- intersect(ecm_genes, rownames(group_means))

cat("ECM genes found in dataset:", length(ecm_present), "\n")
print(ecm_present)

# Rank these by variance across groups and take top 25
ecm_var <- apply(group_means[ecm_present, , drop = FALSE], 1, var)
top25_ecm <- names(sort(ecm_var, decreasing = TRUE))[1:min(25, length(ecm_present))]

top_mat_ecm <- group_means[top25_ecm, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_ecm,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 11,
              cellheight = 14,
              cellwidth = 70,
              treeheight_row = 30,
              border_color = "white",
              main = "Top 25 ECM Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_ecm)
n_cols <- ncol(top_mat_ecm)
pdf_height <- (n_rows * 14) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("ecm_top25_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()



# Curated set of catabolic-process-associated genes
# (covers proteasome, lysosomal proteases, autophagy, and general catabolic/proteolytic machinery)
catabolic_genes <- c(
  # Proteasome subunits
  "PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7",
  "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7",
  "PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6",
  "PSMD1","PSMD2","PSMD3","PSMD4","PSMD11","PSMD12","PSMD14",
  # Ubiquitin-proteasome system
  "UBA1","UBE2D1","UBE2N","UBE3A","UBB","UBC",
  # Lysosomal proteases / cathepsins
  "CTSA","CTSB","CTSD","CTSK","CTSL","CTSS","CTSZ",
  # Autophagy machinery
  "LAMP1","LAMP2","SQSTM1","MAP1LC3B","ATG5","ATG7","ATG12","BECN1",
  # General proteases / peptidases
  "LONP1","CAPN1","CAPN2","MMP2","MMP9","MMP14",
  # Metabolic catabolism (protein/lipid/carb breakdown)
  "HSPA8","HSPA9","GAA","NPC1","ACAA2","ACADVL"
)

# Intersect with genes actually present in your data
catabolic_present <- intersect(catabolic_genes, rownames(group_means))

cat("Catabolic genes found in dataset:", length(catabolic_present), "\n")
print(catabolic_present)

# Rank these by variance across groups and take top 10
cat_var <- apply(group_means[catabolic_present, , drop = FALSE], 1, var)
top10_cat <- names(sort(cat_var, decreasing = TRUE))[1:min(10, length(catabolic_present))]

top_mat_cat <- group_means[top10_cat, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_cat,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 11,
              cellheight = 20,
              cellwidth = 70,
              treeheight_row = 25,
              border_color = "white",
              main = "Top 10 Catabolic Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_cat)
n_cols <- ncol(top_mat_cat)
pdf_height <- (n_rows * 20) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("catabolic_top10_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()


# Curated set of catabolic-process-associated genes
# (proteasome, ubiquitin system, lysosomal proteases, autophagy, general proteolysis/metabolism)
catabolic_genes <- c(
  # Proteasome subunits
  "PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7",
  "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMB10",
  "PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6",
  "PSMD1","PSMD2","PSMD3","PSMD4","PSMD6","PSMD7","PSMD8","PSMD11","PSMD12","PSMD13","PSMD14",
  # Ubiquitin-proteasome system
  "UBA1","UBA52","UBE2D1","UBE2D2","UBE2N","UBE2L3","UBE3A","UBB","UBC",
  # Lysosomal proteases / cathepsins
  "CTSA","CTSB","CTSD","CTSK","CTSL","CTSS","CTSZ","CTSC","CTSF","CTSH",
  # Autophagy machinery
  "LAMP1","LAMP2","SQSTM1","MAP1LC3B","ATG3","ATG5","ATG7","ATG9A","ATG12","BECN1","GABARAP",
  # General proteases / peptidases
  "LONP1","CAPN1","CAPN2","CAPNS1","MMP2","MMP9","MMP14","ADAM17","USP7","USP14",
  # Metabolic catabolism (protein/lipid/carb breakdown)
  "HSPA8","HSPA9","GAA","NPC1","ACAA2","ACADVL","ACOX1","HADHA","HADHB","IDH3A"
)

# Intersect with genes actually present in your data
catabolic_present <- intersect(catabolic_genes, rownames(group_means))

cat("Catabolic genes found in dataset:", length(catabolic_present), "\n")
print(catabolic_present)

# Rank these by variance across groups and take top 25
cat_var <- apply(group_means[catabolic_present, , drop = FALSE], 1, var)
top25_cat <- names(sort(cat_var, decreasing = TRUE))[1:min(25, length(catabolic_present))]

top_mat_cat <- group_means[top25_cat, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_cat,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 11,
              cellheight = 14,
              cellwidth = 70,
              treeheight_row = 30,
              border_color = "white",
              main = "Top 25 Catabolic Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_cat)
n_cols <- ncol(top_mat_cat)
pdf_height <- (n_rows * 14) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("catabolic_top25_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()



# Curated set of fibrosis-associated genes
# (collagen/ECM deposition, myofibroblast markers, pro-fibrotic signaling, TGF-beta pathway)
fibrotic_genes <- c(
  # Collagens (core fibrotic ECM)
  "COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL6A1","COL6A2","COL6A3",
  # Myofibroblast / activation markers
  "ACTA2","TAGLN","VIM","FN1","S100A4",
  # Pro-fibrotic signaling
  "TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","CTGF","PDGFRB",
  # ECM crosslinking / remodeling
  "LOX","LOXL1","LOXL2","TGM2","PLOD1","PLOD2",
  # ECM regulatory / matricellular proteins
  "SPARC","POSTN","THBS1","THBS2","FBN1","FBLN5","BGN","DCN","LUM","VCAN",
  # Fibrosis-associated proteases/inhibitors
  "MMP2","MMP9","TIMP1","TIMP2","SERPINE1","SERPINE2",
  # Inflammatory/fibrotic mediators
  "IL6","TNF","CCL2","PDGFA","PDGFB","VEGFA"
)

# Intersect with genes actually present in your data
fibrotic_present <- intersect(fibrotic_genes, rownames(group_means))

cat("Fibrotic genes found in dataset:", length(fibrotic_present), "\n")
print(fibrotic_present)

# Rank these by variance across groups and take top 10
fib_var <- apply(group_means[fibrotic_present, , drop = FALSE], 1, var)
top10_fib <- names(sort(fib_var, decreasing = TRUE))[1:min(10, length(fibrotic_present))]

top_mat_fib <- group_means[top10_fib, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_fib,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 11,
              cellheight = 20,
              cellwidth = 70,
              treeheight_row = 25,
              border_color = "white",
              main = "Top 10 Fibrotic Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_fib)
n_cols <- ncol(top_mat_fib)
pdf_height <- (n_rows * 20) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("fibrotic_top10_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()


# Curated set of fibrosis-associated genes
# (collagen/ECM deposition, myofibroblast markers, TGF-beta signaling, crosslinking, matricellular proteins)
fibrotic_genes <- c(
  # Collagens (core fibrotic ECM)
  "COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL8A1",
  # Myofibroblast / activation markers
  "ACTA2","TAGLN","VIM","FN1","S100A4","CNN1","MYH11",
  # Pro-fibrotic signaling
  "TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7","CTGF","PDGFRB","PDGFRA",
  # ECM crosslinking / remodeling
  "LOX","LOXL1","LOXL2","LOXL3","TGM2","PLOD1","PLOD2","PLOD3",
  # ECM regulatory / matricellular proteins
  "SPARC","POSTN","THBS1","THBS2","FBN1","FBN2","FBLN5","FBLN1","BGN","DCN","LUM","VCAN",
  # Fibrosis-associated proteases/inhibitors
  "MMP2","MMP9","MMP14","TIMP1","TIMP2","TIMP3","SERPINE1","SERPINE2",
  # Inflammatory/fibrotic mediators
  "IL6","IL13","TNF","CCL2","CCL3","PDGFA","PDGFB","VEGFA","CTSK"
)

# Intersect with genes actually present in your data
fibrotic_present <- intersect(fibrotic_genes, rownames(group_means))

cat("Fibrotic genes found in dataset:", length(fibrotic_present), "\n")
print(fibrotic_present)

# Rank these by variance across groups and take top 25
fib_var <- apply(group_means[fibrotic_present, , drop = FALSE], 1, var)
top25_fib <- names(sort(fib_var, decreasing = TRUE))[1:min(25, length(fibrotic_present))]

top_mat_fib <- group_means[top25_fib, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_fib,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 11,
              cellheight = 14,
              cellwidth = 70,
              treeheight_row = 30,
              border_color = "white",
              main = "Top 25 Fibrotic Proteins (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_fib)
n_cols <- ncol(top_mat_fib)
pdf_height <- (n_rows * 14) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("fibrotic_top25_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()




################ fibrotic gene##########################

# Custom gene list (10 genes)
custom_genes <- c("BGN","TIMP1","COL3A1","PLOD3","COL1A2",
                  "COL6A1","COL6A3","COL6A2","SPARC","TAGLN")

# Intersect with genes actually present in your data
custom_present <- custom_genes[custom_genes %in% rownames(group_means)]

cat("Genes found in dataset:", length(custom_present), "out of", length(custom_genes), "\n")
print(custom_present)

top_mat_custom <- group_means[custom_present, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_custom,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 12,
              cellheight = 18,
              cellwidth = 75,
              treeheight_row = 30,
              border_color = "white",
              main = "Fibrotic/ECM Marker Panel (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 12)

n_rows <- nrow(top_mat_custom)
n_cols <- ncol(top_mat_custom)
pdf_height <- (n_rows * 18) / 72 + 1.5
pdf_width  <- (n_cols * 75) / 72 + 3.5

pdf("custom_panel_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()

###### sNPC vs sNPC-MG###########################

# Custom gene list (10 genes)
custom_genes <- c("BGN","TIMP1","COL3A1","PLOD3","COL1A2",
                  "COL6A1","COL6A3","COL6A2","SPARC","TAGLN")

custom_present <- custom_genes[custom_genes %in% rownames(group_means)]
cat("Genes found in dataset:", length(custom_present), "out of", length(custom_genes), "\n")
print(custom_present)

# Select only sNPC and sNPC-MG columns
groups_to_keep <- c("sNPC", "sNPC-MG")   # adjust names to exactly match your column names
groups_to_keep <- groups_to_keep[groups_to_keep %in% colnames(group_means)]

cat("Groups found:", groups_to_keep, "\n")

top_mat_custom <- group_means[custom_present, groups_to_keep, drop = FALSE]

# Subset annotation to match the selected groups
annotation_col_sub <- annotation_col_grp[groups_to_keep, , drop = FALSE]

# Plot heatmap
p <- pheatmap(top_mat_custom,
              scale = "row",
              annotation_col = annotation_col_sub,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 12,
              cellheight = 18,
              cellwidth = 90,        # wider since only 2 columns now
              treeheight_row = 30,
              border_color = "white",
              main = "Fibrotic/ECM Marker Panel (sNPC vs sNPC-MG)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 12)

n_rows <- nrow(top_mat_custom)
n_cols <- ncol(top_mat_custom)
pdf_height <- (n_rows * 18) / 72 + 1.5
pdf_width  <- (n_cols * 90) / 72 + 3.5

pdf("custom_panel_sNPC_vs_sNPCMG.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()

###############NP phenotype#################

# Curated set of nucleus pulposus (NP) phenotypic marker genes
# (classic NP-specific markers from disc biology literature - notochordal & mature NP markers)
np_genes <- c(
  "ACAN","COL2A1","SOX9","KRT19","KRT8","KRT18",
  "CA3","FOXF1","PAX1","SHH","TBXT","CD24",
  "GDF5","HIF1A","SLC2A1","VCAN","CDH2","MGP",
  "S100A1","ITGA6","COL9A1","COL11A1","SOX5","SOX6"
)

# Intersect with genes actually present in your data
np_present <- intersect(np_genes, rownames(group_means))

cat("NP phenotypic genes found in dataset:", length(np_present), "\n")
print(np_present)

# Rank by variance across groups and take top 10
np_var <- apply(group_means[np_present, , drop = FALSE], 1, var)
top10_np <- names(sort(np_var, decreasing = TRUE))[1:min(10, length(np_present))]

top_mat_np <- group_means[top10_np, , drop = FALSE]

# Plot heatmap across ALL groups
p <- pheatmap(top_mat_np,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 10,
              fontsize_col = 11,
              cellheight = 20,
              cellwidth = 70,
              treeheight_row = 25,
              border_color = "white",
              main = "Top 10 NP Phenotypic Markers (Group Averages)",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 10)
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 11)

n_rows <- nrow(top_mat_np)
n_cols <- ncol(top_mat_np)
pdf_height <- (n_rows * 20) / 72 + 1.5
pdf_width  <- (n_cols * 70) / 72 + 3.5

pdf("NP_phenotype_top10_heatmap.pdf", width = pdf_width, height = pdf_height)
grid.newpage()
grid.draw(p$gtable)
dev.off()
####################sNPC vs sNPC-MG catabolic marker##################

# Curated catabolic protein panel (matrix-degrading enzymes / catabolic pathway proteins)
catabolic_proteins <- c(
  "MMP1","MMP2","MMP3","MMP7","MMP8","MMP9","MMP10","MMP12","MMP13","MMP14",
  "ADAMTS1","ADAMTS4","ADAMTS5","ADAMTS9",
  "CTSK","CTSB","CTSD","CTSL","CTSS",
  "PLAU","PLAT","PLAUR",
  "TIMP1","TIMP2","TIMP3","TIMP4"
)

# --- Robust case-insensitive matching against your dataset ---
rn <- rownames(group_means)
rn_upper <- toupper(trimws(rn))
proteins_upper <- toupper(catabolic_proteins)

matched_idx <- match(proteins_upper, rn_upper)
found <- !is.na(matched_idx)

matched_proteins_original_name <- rn[matched_idx[found]]
matched_query_name <- catabolic_proteins[found]

cat("Proteins found:", length(matched_proteins_original_name), "out of", length(catabolic_proteins), "\n")
print(data.frame(query = matched_query_name, matched_rowname = matched_proteins_original_name))

matched_proteins_unique <- unique(matched_proteins_original_name)

# --- Select only sNPC and sNPC-MG columns ---
groups_to_keep <- c("sNPC", "sNPC-MG")   # adjust to exactly match your column names
groups_to_keep <- groups_to_keep[groups_to_keep %in% colnames(group_means)]
cat("Groups found:", groups_to_keep, "\n")

if (length(groups_to_keep) < 2) {
  cat("WARNING: Less than 2 matching groups found. Check colnames(group_means):\n")
  print(colnames(group_means))
}

top_mat <- group_means[matched_proteins_unique, groups_to_keep, drop = FALSE]

# --- Drop NA / zero-variance rows (prevents scaling errors) ---
row_sd <- apply(top_mat, 1, sd, na.rm = TRUE)
valid_proteins <- names(row_sd)[!is.na(row_sd) & row_sd > 0]
top_mat_valid <- top_mat[valid_proteins, , drop = FALSE]

cat("Proteins with valid variance:", nrow(top_mat_valid), "\n")
print(rownames(top_mat_valid))

# --- Take top 25 by variance (or fewer if less available) ---
protein_var <- apply(top_mat_valid, 1, var)
top25 <- names(sort(protein_var, decreasing = TRUE))[1:min(25, length(protein_var))]
top_mat_final <- top_mat_valid[top25, , drop = FALSE]

# --- Subset annotation to match selected groups ---
annotation_col_sub <- annotation_col_grp[groups_to_keep, , drop = FALSE]

# --- Plot heatmap (only if at least 2 rows available) ---
if (nrow(top_mat_final) >= 2) {
  p <- pheatmap(top_mat_final,
                scale = "row",
                annotation_col = annotation_col_sub,
                annotation_colors = ann_colors,
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                fontsize_row = 9,
                fontsize_col = 12,
                cellheight = 14,
                cellwidth = 90,
                treeheight_row = 30,
                border_color = "white",
                main = "Top 25 Catabolic Proteins (sNPC vs sNPC-MG)",
                silent = TRUE)
  
  p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp <- gpar(fontface = "bold", fontsize = 9)
  p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "bold", fontsize = 12)
  
  n_rows <- nrow(top_mat_final)
  n_cols <- ncol(top_mat_final)
  pdf_height <- (n_rows * 14) / 72 + 1.5
  pdf_width  <- (n_cols * 90) / 72 + 3.5
  
  pdf("catabolic_protein_top25_sNPC_vs_sNPCMG.pdf", width = pdf_width, height = pdf_height)
  grid.newpage()
  grid.draw(p$gtable)
  dev.off()
  
  cat("Heatmap saved successfully.\n")
} else {
  cat("ERROR: Fewer than 2 valid proteins found. Check matching output above.\n")
}


######################heat map for mitochondrial protein###############

# Mitochondrial stress / UPRmt (mitochondrial unfolded protein response) protein panel
mito_stress_genes <- c(
  "HSPD1","HSPE1","HSPA9","LONP1","CLPP","CLPX","YME1L1","AFG3L2","OMA1",
  "OPA1","DNM1L","MFN1","MFN2","PINK1","PRKN","SOD2","UCP2","TFAM","POLG",
  "ATF4","ATF5","DDIT3","HSPA1A","DNAJA3","HSPB1"
)

# --- Robust case-insensitive matching against group_means ---
rn <- rownames(group_means)
rn_upper <- toupper(trimws(rn))
genes_upper <- toupper(mito_stress_genes)

matched_idx <- match(genes_upper, rn_upper)
found <- !is.na(matched_idx)

matched_original_name <- rn[matched_idx[found]]
matched_query_name <- mito_stress_genes[found]

cat("Genes found:", length(matched_original_name), "out of", length(mito_stress_genes), "\n")
print(data.frame(query = matched_query_name, matched_rowname = matched_original_name))

matched_unique <- unique(matched_original_name)
top_mat_grp <- group_means[matched_unique, , drop = FALSE]

# --- Drop NA / zero-variance rows ---
row_sd <- apply(top_mat_grp, 1, sd, na.rm = TRUE)
valid_genes <- names(row_sd)[!is.na(row_sd) & row_sd > 0]
top_mat_valid <- top_mat_grp[valid_genes, , drop = FALSE]

cat("Genes with valid variance:", nrow(top_mat_valid), "\n")
print(rownames(top_mat_valid))

# --- Plot heatmap ---
pheatmap(top_mat_valid,
         scale = "row",
         annotation_col = annotation_col_grp,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 9,
         fontsize_col = 11,
         cellheight = 14,
         cellwidth = 70,
         treeheight_row = 30,
         border_color = "white",
         main = "Mitochondrial Stress Proteins (Group Averages)",
         filename = "mito_stress_group_averages_heatmap.pdf",
         width = ncol(top_mat_valid) * 1.0 + 3.5,
         height = nrow(top_mat_valid) * 0.2 + 1.5)

cat("Heatmap saved successfully.\n")]

########### Mitochondria and oxphos###################


# ============================================================
# Comprehensive Mitochondria + OXPHOS protein panel
# ============================================================

oxphos_mito_genes <- c(
  # Complex I (NADH:ubiquinone oxidoreductase)
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
  "NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFAB1",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFV3",
  
  # Complex II (Succinate dehydrogenase)
  "SDHA","SDHB","SDHC","SDHD",
  
  # Complex III (Cytochrome bc1 complex)
  "UQCRC1","UQCRC2","UQCRB","UQCRQ","UQCRH","UQCRFS1","CYC1","CYCS",
  
  # Complex IV (Cytochrome c oxidase)
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7B","COX7C","COX8A",
  "MT-CO1","MT-CO2","MT-CO3",
  
  # Complex V (ATP synthase)
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5PB","ATP5PD","ATP5PO","ATP5MC1","ATP5ME","ATP5MF","ATP5MG",
  
  # Mitochondrial structure / transport
  "TOMM20","TOMM22","TIMM23","TIMM50","VDAC1","VDAC2","VDAC3",
  
  # TCA cycle / metabolism
  "CS","IDH2","IDH3A","MDH2","ACO2","FH","SUCLA2","PDHA1","PDHB",
  
  # Mitochondrial biogenesis / dynamics / stress
  "TFAM","POLG","PPARGC1A","OPA1","MFN1","MFN2","DNM1L",
  "SOD2","PINK1","HSPD1","HSPE1"
)

# --- Robust case-insensitive matching against group_means ---
rn <- rownames(group_means)
rn_upper <- toupper(trimws(rn))
genes_upper <- toupper(oxphos_mito_genes)

matched_idx <- match(genes_upper, rn_upper)
found <- !is.na(matched_idx)

matched_original_name <- rn[matched_idx[found]]
matched_query_name <- oxphos_mito_genes[found]

cat("Genes found:", length(matched_original_name), "out of", length(oxphos_mito_genes), "\n")
print(data.frame(query = matched_query_name, matched_rowname = matched_original_name))

matched_unique <- unique(matched_original_name)
top_mat <- group_means[matched_unique, , drop = FALSE]

# --- Drop NA / zero-variance rows ---
row_sd <- apply(top_mat, 1, sd, na.rm = TRUE)
valid_genes <- names(row_sd)[!is.na(row_sd) & row_sd > 0]
top_mat_valid <- top_mat[valid_genes, , drop = FALSE]

cat("Genes with valid variance:", nrow(top_mat_valid), "\n")

# --- If too many rows, keep top 40 by variance for readability ---
max_rows <- 40
if (nrow(top_mat_valid) > max_rows) {
  gene_var <- apply(top_mat_valid, 1, var)
  keep <- names(sort(gene_var, decreasing = TRUE))[1:max_rows]
  top_mat_final <- top_mat_valid[keep, , drop = FALSE]
  cat("Trimmed to top", max_rows, "by variance for plotting.\n")
} else {
  top_mat_final <- top_mat_valid
}

print(rownames(top_mat_final))

# --- Plot heatmap ---
pheatmap(top_mat_final,
         scale = "row",
         annotation_col = annotation_col_grp,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 8,
         fontsize_col = 11,
         cellheight = 12,
         cellwidth = 70,
         treeheight_row = 30,
         border_color = "white",
         main = "Mitochondrial / OXPHOS Proteins (Group Averages)",
         filename = "mito_oxphos_group_averages_heatmap.pdf",
         width = ncol(top_mat_final) * 1.0 + 3.5,
         height = nrow(top_mat_final) * 0.18 + 1.5)

cat("Heatmap saved successfully.\n")

# --- Step 1: Create heatmap WITHOUT fontface argument inside pheatmap() ---
p <- pheatmap(top_mat_final,
              scale = "row",
              annotation_col = annotation_col_grp,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 11,
              cellheight = 14,
              cellwidth = 70,
              border_color = "white",
              main = "Mitochondrial / OXPHOS Proteins (Group Averages)",
              silent = TRUE)

# --- Step 2: Force bold AFTER creation (this is the only reliable way) ---
idx_row <- which(p$gtable$layout$name == "row_names")
if (length(idx_row) > 0) {
  p$gtable$grobs[[idx_row]]$gp <- grid::gpar(fontface = "bold", fontsize = 9)
}

idx_col <- which(p$gtable$layout$name == "col_names")
if (length(idx_col) > 0) {
  p$gtable$grobs[[idx_col]]$gp <- grid::gpar(fontface = "bold", fontsize = 11)
}

idx_title <- which(p$gtable$layout$name == "main")
if (length(idx_title) > 0) {
  p$gtable$grobs[[idx_title]]$gp <- grid::gpar(fontface = "bold", fontsize = 13)
}

# --- Step 3: Save as PDF ---
pdf("heatmap_bold.pdf", width = 8, height = 6)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

cat("Bold heatmap saved successfully.\n")

################################sNPC vs sNPC-MG mitochindria ##############

# ============================================================
# Mitochondria + OXPHOS protein panel — sNPC vs sNPC-MG only (bold text)
# ============================================================

oxphos_mito_genes <- c(
  # Complex I
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
  "NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFAB1",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFV3",
  
  # Complex II
  "SDHA","SDHB","SDHC","SDHD",
  
  # Complex III
  "UQCRC1","UQCRC2","UQCRB","UQCRQ","UQCRH","UQCRFS1","CYC1","CYCS",
  
  # Complex IV
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7B","COX7C","COX8A",
  "MT-CO1","MT-CO2","MT-CO3",
  
  # Complex V
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5PB","ATP5PD","ATP5PO","ATP5MC1","ATP5ME","ATP5MF","ATP5MG",
  
  # Mitochondrial structure/transport
  "TOMM20","TOMM22","TIMM23","TIMM50","VDAC1","VDAC2","VDAC3",
  
  # TCA cycle / metabolism
  "CS","IDH2","IDH3A","MDH2","ACO2","FH","SUCLA2","PDHA1","PDHB",
  
  # Biogenesis / dynamics / stress
  "TFAM","POLG","PPARGC1A","OPA1","MFN1","MFN2","DNM1L",
  "SOD2","PINK1","HSPD1","HSPE1"
)

# --- Robust case-insensitive matching against group_means ---
rn <- rownames(group_means)
rn_upper <- toupper(trimws(rn))
genes_upper <- toupper(oxphos_mito_genes)

matched_idx <- match(genes_upper, rn_upper)
found <- !is.na(matched_idx)

matched_original_name <- rn[matched_idx[found]]
matched_query_name <- oxphos_mito_genes[found]

cat("Genes found:", length(matched_original_name), "out of", length(oxphos_mito_genes), "\n")
print(data.frame(query = matched_query_name, matched_rowname = matched_original_name))

matched_unique <- unique(matched_original_name)

# --- Select only sNPC and sNPC-MG columns ---
groups_to_keep <- c("sNPC", "sNPC-MG")   # adjust if your colnames differ (e.g. "sNPC_MG")
groups_to_keep <- groups_to_keep[groups_to_keep %in% colnames(group_means)]
cat("Groups found:", groups_to_keep, "\n")

if (length(groups_to_keep) < 2) {
  cat("WARNING: Less than 2 matching groups found. Check colnames(group_means):\n")
  print(colnames(group_means))
}

top_mat <- group_means[matched_unique, groups_to_keep, drop = FALSE]

# --- Drop NA / zero-variance rows ---
row_sd <- apply(top_mat, 1, sd, na.rm = TRUE)
valid_genes <- names(row_sd)[!is.na(row_sd) & row_sd > 0]
top_mat_valid <- top_mat[valid_genes, , drop = FALSE]

cat("Genes with valid variance:", nrow(top_mat_valid), "\n")

# --- If too many rows, keep top 40 by variance for readability ---
max_rows <- 40
if (nrow(top_mat_valid) > max_rows) {
  gene_var <- apply(top_mat_valid, 1, var)
  keep <- names(sort(gene_var, decreasing = TRUE))[1:max_rows]
  top_mat_final <- top_mat_valid[keep, , drop = FALSE]
  cat("Trimmed to top", max_rows, "by variance for plotting.\n")
} else {
  top_mat_final <- top_mat_valid
}

print(rownames(top_mat_final))

# --- Subset annotation to match selected groups ---
annotation_col_sub <- annotation_col_grp[groups_to_keep, , drop = FALSE]

# --- Step 1: Create heatmap (NO fontface argument inside pheatmap) ---
p <- pheatmap(top_mat_final,
              scale = "row",
              annotation_col = annotation_col_sub,
              annotation_colors = ann_colors,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              fontsize_row = 9,
              fontsize_col = 12,
              cellheight = 14,
              cellwidth = 90,
              treeheight_row = 30,
              border_color = "white",
              main = "Mitochondrial / OXPHOS Proteins (sNPC vs sNPC-MG)",
              silent = TRUE)

# --- Step 2: Force bold on row names, column names, and title ---
idx_row <- which(p$gtable$layout$name == "row_names")
if (length(idx_row) > 0) {
  p$gtable$grobs[[idx_row]]$gp <- grid::gpar(fontface = "bold", fontsize = 9)
}

idx_col <- which(p$gtable$layout$name == "col_names")
if (length(idx_col) > 0) {
  p$gtable$grobs[[idx_col]]$gp <- grid::gpar(fontface = "bold", fontsize = 12)
}

idx_title <- which(p$gtable$layout$name == "main")
if (length(idx_title) > 0) {
  p$gtable$grobs[[idx_title]]$gp <- grid::gpar(fontface = "bold", fontsize = 13)
}

# --- Step 3: Save as PDF ---
n_rows <- nrow(top_mat_final)
n_cols <- ncol(top_mat_final)
pdf_height <- (n_rows * 14) / 72 + 1.5
pdf_width  <- (n_cols * 90) / 72 + 3.5

pdf("mito_oxphos_sNPC_vs_sNPCMG_bold.pdf", width = pdf_width, height = pdf_height)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

cat("Bold heatmap saved successfully.\n")

##################### venn diagram#########


# ============================================================
# Venn Diagram — compare protein sets across groups
# ============================================================
install.packages("VennDiagram")


library(VennDiagram)
library(grid)

# --- Define protein sets per group ---
# Example: proteins detected (non-NA, above threshold) in each group
# Adjust the criteria below to match what defines "present" in your data
# (e.g., detected/non-missing, or significantly changed vs control)

get_detected_proteins <- function(mat, group_col, threshold = 0) {
  vals <- mat[, group_col]
  rownames(mat)[!is.na(vals) & vals > threshold]
}

# --- Example for 2 groups: sNPC vs sNPC-MG ---
set_sNPC    <- get_detected_proteins(group_means, "sNPC")
set_sNPC_MG <- get_detected_proteins(group_means, "sNPC-MG")

venn.plot <- venn.diagram(
  x = list(sNPC = set_sNPC, `sNPC-MG` = set_sNPC_MG),
  filename = NULL,
  fill = c("#4E79A7", "#F28E2B"),
  alpha = 0.6,
  cex = 1.3,
  cat.cex = 1.3,
  cat.fontface = "bold",
  fontface = "bold",
  main = "Protein Overlap: sNPC vs sNPC-MG",
  main.fontface = "bold",
  main.cex = 1.4,
  lwd = 2,
  col = "white"
)

pdf("venn_sNPC_vs_sNPC-MG.pdf", width = 6, height = 6)
grid.newpage()
grid.draw(venn.plot)
dev.off()

cat("Venn diagram saved successfully.\n")



########################### protein distirbution across sample#########


# ============================================================
# Protein Distribution Across Groups and Samples
# ============================================================
library(ggplot2)
library(reshape2)   # for melting data into long format
library(RColorBrewer)

# --- Assumptions ---
# 'data_mat'  : matrix/data.frame of proteins (rows) x samples (columns), values = intensity/abundance
# 'sample_metadata' : data.frame with columns "sample" and "group", mapping each sample to its group
#
# Adjust variable names below if yours differ (e.g., raw intensity matrix, expr_matrix, etc.)

# --- Step 1: Convert wide matrix to long format ---
data_long <- melt(as.matrix(data_mat), varnames = c("Protein", "Sample"), value.name = "Intensity")

# --- Step 2: Merge with group info ---
data_long <- merge(data_long, sample_metadata, by.x = "Sample", by.y = "sample")

# --- Step 3: Remove NA/zero values (optional, avoids skewed plots) ---
data_long <- data_long[!is.na(data_long$Intensity) & data_long$Intensity > 0, ]

# --- Step 4: Log-transform if data is raw intensity (skip if already log-scale) ---
data_long$log_Intensity <- log2(data_long$Intensity + 1)

# --- Step 5: Order samples by group for cleaner x-axis ---
data_long$Sample <- factor(data_long$Sample,
                           levels = sample_metadata$sample[order(sample_metadata$group)])

# --- Step 6: Plot — Boxplot per sample, colored by group ---
n_groups <- length(unique(data_long$group))
group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]

p <- ggplot(data_long, aes(x = Sample, y = log_Intensity, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Protein Intensity Distribution Across Samples",
       x = "Sample", y = "log2(Intensity)", fill = "Group") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

print(p)

ggsave("protein_distribution_samples.pdf", plot = p, width = 10, height = 6)
cat("Plot saved as protein_distribution_samples.pdf\n")

# ============================================================
# Optional: Distribution summarized by GROUP only (not per sample)
# ============================================================
p_group <- ggplot(data_long, aes(x = group, y = log_Intensity, fill = group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Protein Intensity Distribution by Group",
       x = "Group", y = "log2(Intensity)") +
  theme_bw(base_size = 13) +
  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

print(p_group)
ggsave("protein_distribution_groups.pdf", plot = p_group, width = 7, height = 6)
cat("Plot saved as protein_distribution_groups.pdf\n")


################### specific protein#################


# ============================================================
# Specific Protein Panel — Bar graphs with Mean ± SD across groups
# Proteins: IL1B, MMP3, IL8 (CXCL8), HSPA5/GRP78, DDIT3/CHOP, ASIC1, ASIC2, ASIC3
# ============================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# --- Target protein panel (using standard gene symbols) ---
target_genes <- c(
  "IL1B",    # IL-1 beta
  "MMP3",    # MMP-3
  "CXCL8",   # IL-8 (alt symbol: IL8)
  "IL8",     # backup alt symbol
  "HSPA5",   # GRP78 / BiP (ER stress) -- assumed match for "CGRP78"
  "DDIT3",   # CHOP
  "ASIC1",
  "ASIC2",
  "ASIC3"
)

# --- Match against your data (case-insensitive) ---
rn <- rownames(data_mat)          # <-- REPLACE with your actual protein x sample matrix name
rn_upper <- toupper(trimws(rn))
genes_upper <- toupper(target_genes)

matched_idx <- match(genes_upper, rn_upper)
found <- !is.na(matched_idx)

matched_rows <- rn[matched_idx[found]]
matched_query <- target_genes[found]

cat("Proteins found in data:\n")
print(data.frame(query = matched_query, matched_rowname = matched_rows))

if (length(matched_rows) == 0) {
  stop("No matching proteins found. Check gene symbols against rownames(data_mat).")
}

matched_unique <- unique(matched_rows)

# --- Subset data to just these proteins ---
sub_mat <- data_mat[matched_unique, , drop = FALSE]   # <-- REPLACE data_mat if needed

# --- Convert to long format ---
data_long <- melt(as.matrix(sub_mat), varnames = c("Protein", "Sample"), value.name = "Intensity")

# --- Merge with sample -> group mapping ---
data_long <- merge(data_long, sample_metadata, by.x = "Sample", by.y = "sample")  # <-- REPLACE sample_metadata

# --- Remove missing values ---
data_long <- data_long[!is.na(data_long$Intensity), ]

# --- Compute mean, SD, SEM per Protein x Group ---
summary_df <- data_long %>%
  group_by(Protein, group) %>%
  summarise(
    Mean = mean(Intensity, na.rm = TRUE),
    SD   = sd(Intensity, na.rm = TRUE),
    N    = sum(!is.na(Intensity)),
    SEM  = SD / sqrt(N),
    .groups = "drop"
  )

print(summary_df)

# --- Plot: one bar chart per protein (faceted), mean +/- SD ---
n_groups <- length(unique(summary_df$group))
group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]

p <- ggplot(summary_df, aes(x = group, y = Mean, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.25, linewidth = 0.7) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Protein Expression Across Groups (Mean ± SD)",
       x = "Group", y = "Intensity") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

print(p)

ggsave("protein_panel_barplot_SD.pdf", plot = p, width = 10, height = 8)
cat("Bar plot saved as protein_panel_barplot_SD.pdf\n")

ls()

# ============================================================
# Auto-detect the correct protein x sample matrix
# ============================================================

candidates <- c("log_mat_imp", "log_mat_filt", "log_mat", "normalized_matrix",
                "expr_mat", "abundance_matrix", "raw_data", "expr_agg")

candidates <- candidates[sapply(candidates, exists)]

cat("Candidate objects found in environment:\n")
for (obj_name in candidates) {
  obj <- get(obj_name)
  if (is.matrix(obj) || is.data.frame(obj)) {
    cat(sprintf(" - %-20s dim: %d rows x %d cols\n", obj_name, nrow(obj), ncol(obj)))
  }
}

# Check sample_metadata structure
cat("\nsample_metadata preview:\n")
print(head(sample_metadata))
cat("Columns in sample_metadata:", paste(colnames(sample_metadata), collapse=", "), "\n")

for (obj_name in c("log_mat_imp","log_mat_filt","log_mat","normalized_matrix","expr_mat","abundance_matrix","raw_data","expr_agg")) {
  if (exists(obj_name)) {
    obj <- get(obj_name)
    cat(obj_name, ":", nrow(obj), "rows x", ncol(obj), "cols | sample cols:", paste(head(colnames(obj)), collapse=", "), "\n")
  }
}


# ============================================================
# Specific Protein Panel — Bar graphs with Mean ± SD across groups
# Using: log_mat_imp (filtered, imputed matrix)
# ============================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# --- Target protein panel ---
target_genes <- c(
  "IL1B", "MMP3", "CXCL8", "IL8",   # IL8 alt symbol as backup
  "HSPA5",   # GRP78/BiP -- assumed match for "CGRP78"
  "DDIT3",   # CHOP
  "ASIC1", "ASIC2", "ASIC3"
)

# --- Match against log_mat_imp rownames (case-insensitive) ---
rn <- rownames(log_mat_imp)
rn_upper <- toupper(trimws(rn))
genes_upper <- toupper(target_genes)

matched_idx <- match(genes_upper, rn_upper)
found <- !is.na(matched_idx)

matched_rows <- rn[matched_idx[found]]
matched_query <- target_genes[found]

cat("Proteins found in data:\n")
print(data.frame(query = matched_query, matched_rowname = matched_rows))

missing_genes <- target_genes[!found]
if (length(missing_genes) > 0) {
  cat("\nNOT found in log_mat_imp rownames:", paste(missing_genes, collapse = ", "), "\n")
}

if (length(matched_rows) == 0) stop("No matching proteins found. Check gene symbols.")

matched_unique <- unique(matched_rows)

# --- Subset matrix to just these proteins ---
sub_mat <- log_mat_imp[matched_unique, , drop = FALSE]

# --- Convert to long format ---
data_long <- melt(as.matrix(sub_mat), varnames = c("Protein", "Sample"), value.name = "log2_Intensity")

# --- Merge with sample -> group mapping ---
data_long <- merge(data_long, sample_metadata[, c("Sample", "Group")], by = "Sample")

# --- Remove missing values ---
data_long <- data_long[!is.na(data_long$log2_Intensity), ]

# --- Compute mean, SD, SEM per Protein x Group ---
summary_df <- data_long %>%
  group_by(Protein, Group) %>%
  summarise(
    Mean = mean(log2_Intensity, na.rm = TRUE),
    SD   = sd(log2_Intensity, na.rm = TRUE),
    N    = sum(!is.na(log2_Intensity)),
    SEM  = SD / sqrt(N),
    .groups = "drop"
  )

print(summary_df)

# --- Plot: bar chart per protein (faceted), mean +/- SD ---
n_groups <- length(unique(summary_df$Group))
group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]

p <- ggplot(summary_df, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.25, linewidth = 0.7) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Protein Expression Across Groups (Mean ± SD)",
       x = "Group", y = "log2 Intensity") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

print(p)

ggsave("protein_panel_barplot_SD.pdf", plot = p, width = 10, height = 8)
cat("\nBar plot saved as protein_panel_barplot_SD.pdf\n")



################## KEGG pathway sNPC vs sNPC-MGt##############

# ============================================================
# STEP 1: Required packages
# ============================================================
required_pkgs <- c("limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ============================================================
# STEP 2: Check group levels in metadata
# ============================================================
print(table(sample_metadata$Group))



# ============================================================
# STEP 3: Subset data to sNPC vs sNPC-MG samples
# ============================================================
keep_samples <- sample_metadata$Sample[sample_metadata$Group %in% c("sNPC", "sNPC-MG")]
keep_samples <- intersect(keep_samples, colnames(log_mat_imp))

expr_sub <- log_mat_imp[, keep_samples, drop = FALSE]
meta_sub <- sample_metadata[sample_metadata$Sample %in% keep_samples, ]
meta_sub <- meta_sub[match(colnames(expr_sub), meta_sub$Sample), ]  # align order

meta_sub$Group <- factor(meta_sub$Group, levels = c("sNPC", "sNPC-MG"))

cat("Samples used:\n")
print(data.frame(Sample = meta_sub$Sample, Group = meta_sub$Group))

# ============================================================
# STEP 4: Differential expression with limma (moderated t-test)
# ============================================================
design <- model.matrix(~ Group, data = meta_sub)
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

# Coefficient = sNPC-MG vs sNPC (since sNPC is reference level)
de_results <- topTable(fit, coef = "GroupsNPC-MG", number = Inf, sort.by = "P")
de_results$Protein <- rownames(de_results)

cat("\nTop 10 DE results:\n")
print(head(de_results, 10))

# ============================================================
# STEP 5: Filter significant DE genes
# ============================================================
sig_genes <- de_results[de_results$adj.P.Val < 0.05 & abs(de_results$logFC) > 1, ]
cat("\nNumber of significant DE proteins (adj.P < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")

if (nrow(sig_genes) < 5) {
  cat("Few/no significant genes at strict cutoff -- consider relaxing to p < 0.05 (unadjusted) or |log2FC| > 0.5\n")
}

# ============================================================
# STEP 6: Map gene symbols to Entrez IDs
# ============================================================
gene_symbols <- sig_genes$Protein

entrez_map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("\nMapped", nrow(entrez_map), "of", length(gene_symbols), "gene symbols to Entrez IDs\n")

unmapped <- setdiff(gene_symbols, entrez_map$SYMBOL)
if (length(unmapped) > 0) {
  cat("Unmapped symbols:", paste(unmapped, collapse = ", "), "\n")
}

# ============================================================
# STEP 7: KEGG pathway enrichment
# ============================================================
kegg_result <- enrichKEGG(
  gene         = entrez_map$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_df <- as.data.frame(kegg_result)
cat("\nKEGG pathways enriched:", nrow(kegg_df), "\n")
print(head(kegg_df, 15))

write.csv(kegg_df, "KEGG_sNPC_vs_sNPCMG.csv", row.names = FALSE)

# ============================================================
# STEP 8: Plot results
# ============================================================
if (nrow(kegg_df) > 0) {
  p_dot <- dotplot(kegg_result, showCategory = 15, title = "KEGG Pathway Enrichment: sNPC vs sNPC-MG") +
    theme(axis.text.y = element_text(size = 9))
  print(p_dot)
  ggsave("KEGG_dotplot_sNPC_vs_sNPCMG.pdf", plot = p_dot, width = 9, height = 7)
  
  p_bar <- barplot(kegg_result, showCategory = 15, title = "KEGG Pathway Enrichment: sNPC vs sNPC-MG")
  print(p_bar)
  ggsave("KEGG_barplot_sNPC_vs_sNPCMG.pdf", plot = p_bar, width = 9, height = 7)
  
  cat("\nPlots saved: KEGG_dotplot_sNPC_vs_sNPCMG.pdf, KEGG_barplot_sNPC_vs_sNPCMG.pdf\n")
} else {
  cat("\nNo enriched KEGG pathways found at current thresholds.\n")
}

################################IVD repated patway###################


# ============================================================
# IVD-Relevant KEGG Pathway Filtering
# (Run AFTER kegg_result / kegg_df from previous step)
# ============================================================

# --- Curated list of KEGG pathways commonly implicated in IVD degeneration ---
ivd_relevant_pathways <- c(
  "hsa04512",  # ECM-receptor interaction
  "hsa04510",  # Focal adhesion
  "hsa04668",  # TNF signaling pathway
  "hsa04064",  # NF-kappa B signaling pathway
  "hsa04657",  # IL-17 signaling pathway
  "hsa04010",  # MAPK signaling pathway
  "hsa04151",  # PI3K-Akt signaling pathway
  "hsa04060",  # Cytokine-cytokine receptor interaction
  "hsa04210",  # Apoptosis
  "hsa04115",  # p53 signaling pathway
  "hsa04218",  # Cellular senescence
  "hsa04066",  # HIF-1 signaling pathway
  "hsa04310",  # Wnt signaling pathway
  "hsa04350",  # TGF-beta signaling pathway
  "hsa04141",  # Protein processing in endoplasmic reticulum (ER stress - GRP78/CHOP)
  "hsa05323",  # Rheumatoid arthritis (shared inflammatory ECM breakdown biology)
  "hsa04380",  # Osteoclast differentiation
  "hsa04610",  # Complement and coagulation cascades
  "hsa04611",  # Platelet activation
  "hsa04620",  # Toll-like receptor signaling pathway
  "hsa04621"   # NOD-like receptor signaling pathway
)

# --- Subset your existing KEGG results to just these pathways ---
kegg_ivd_df <- kegg_df[kegg_df$ID %in% ivd_relevant_pathways, ]
kegg_ivd_df <- kegg_ivd_df[order(kegg_ivd_df$p.adjust), ]

cat("IVD-relevant KEGG pathways enriched in sNPC vs sNPC-MG:\n")
print(kegg_ivd_df[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")])

write.csv(kegg_ivd_df, "KEGG_IVD_relevant_sNPC_vs_sNPCMG.csv", row.names = FALSE)

# ============================================================
# Plot: IVD-relevant pathways only
# ============================================================
library(ggplot2)

if (nrow(kegg_ivd_df) > 0) {
  kegg_ivd_df$Description <- factor(kegg_ivd_df$Description,
                                    levels = kegg_ivd_df$Description[order(kegg_ivd_df$p.adjust, decreasing = TRUE)])
  
  p_ivd <- ggplot(kegg_ivd_df, aes(x = Description, y = -log10(p.adjust), size = Count, color = p.adjust)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "IVD-Relevant KEGG Pathways: sNPC vs sNPC-MG",
         x = "", y = "-log10(adjusted p-value)", size = "Gene Count", color = "adj. p-value") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  print(p_ivd)
  ggsave("KEGG_IVD_relevant_dotplot.pdf", plot = p_ivd, width = 9, height = 6)
  cat("\nPlot saved: KEGG_IVD_relevant_dotplot.pdf\n")
} else {
  cat("\nNone of the IVD-relevant pathways were significantly enriched in this comparison.\n")
  cat("Consider: (1) relaxing DE thresholds in Step 5, or (2) checking full kegg_df for near-significant IVD-related hits.\n")
}

# --- Also show full list sorted, highlighting IVD-relevant rows ---
kegg_df$IVD_relevant <- kegg_df$ID %in% ivd_relevant_pathways
cat("\nFull KEGG table (IVD-relevant pathways flagged):\n")
print(kegg_df[order(kegg_df$p.adjust), c("ID", "Description", "p.adjust", "Count", "IVD_relevant")])



# ============================================================
# Top 10 IVD-relevant KEGG pathways - directionality by group
# ============================================================
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

# --- Take top 10 IVD-relevant pathways by significance ---
top10_ivd <- kegg_ivd_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

cat("Top 10 IVD-relevant pathways:\n")
print(top10_ivd[, c("ID", "Description", "p.adjust", "Count")])

# ============================================================
# For each pathway, compute mean log2FC of its DE genes
# (positive logFC = up in sNPC-MG, negative = up in sNPC)
# ============================================================

# Map Entrez IDs in geneID column back to gene symbols
entrez_to_symbol <- setNames(entrez_map$SYMBOL, entrez_map$ENTREZID)

pathway_direction <- lapply(seq_len(nrow(top10_ivd)), function(i) {
  entrez_ids <- strsplit(top10_ivd$geneID[i], "/")[[1]]
  symbols <- entrez_to_symbol[entrez_ids]
  symbols <- symbols[!is.na(symbols)]
  
  fc_values <- sig_genes$logFC[match(symbols, sig_genes$Protein)]
  fc_values <- fc_values[!is.na(fc_values)]
  
  data.frame(
    ID = top10_ivd$ID[i],
    Description = top10_ivd$Description[i],
    p.adjust = top10_ivd$p.adjust[i],
    Count = top10_ivd$Count[i],
    Mean_log2FC = mean(fc_values, na.rm = TRUE),
    Direction = ifelse(mean(fc_values, na.rm = TRUE) > 0, "Up in sNPC-MG", "Up in sNPC")
  )
})

pathway_direction_df <- do.call(rbind, pathway_direction)
pathway_direction_df <- pathway_direction_df[order(pathway_direction_df$p.adjust), ]

cat("\nPathway directionality:\n")
print(pathway_direction_df)

# ============================================================
# Bar graph: -log10(adj.p) with bars colored by group direction
# ============================================================
pathway_direction_df$Description <- factor(
  pathway_direction_df$Description,
  levels = pathway_direction_df$Description[order(pathway_direction_df$p.adjust, decreasing = TRUE)]
)

p_top10_ivd <- ggplot(pathway_direction_df, aes(x = Description, y = -log10(p.adjust), fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Up in sNPC" = "#4575b4", "Up in sNPC-MG" = "#d73027")) +
  labs(title = "Top 10 IVD-Relevant KEGG Pathways: sNPC vs sNPC-MG",
       x = "", y = "-log10(adjusted p-value)", fill = "Direction") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

print(p_top10_ivd)

ggsave("Top10_IVD_KEGG_barplot_byGroup.pdf", plot = p_top10_ivd, width = 9, height = 6)
cat("\nBar plot saved as Top10_IVD_KEGG_barplot_byGroup.pdf\n")


# ============================================================
# DIAGNOSTIC: check how many pathways/genes you actually have
# ============================================================
cat("Total significant DE genes (sig_genes):", nrow(sig_genes), "\n")
cat("Total pathways in kegg_df (from enrichKEGG):", nrow(kegg_df), "\n")
cat("Total pathways overlapping IVD list:", sum(kegg_df$ID %in% ivd_relevant_pathways), "\n\n")

print(kegg_df[, c("ID", "Description", "p.adjust", "Count")])


# ============================================================
# GSEA-based KEGG (uses ALL genes ranked by log2FC, not just DE-significant ones)
# ============================================================
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Map ALL tested genes (not just significant) to Entrez ---
all_genes_df <- de_results  # full limma results, all proteins tested
all_genes_df$Protein <- rownames(all_genes_df)

entrez_all <- bitr(all_genes_df$Protein, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
all_genes_df <- merge(all_genes_df, entrez_all, by.x = "Protein", by.y = "SYMBOL")

# --- Build ranked gene list (named vector: log2FC, named by Entrez ID) ---
gene_list <- all_genes_df$logFC
names(gene_list) <- all_genes_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

cat("Ranked gene list length:", length(gene_list), "\n")

# --- Run GSEA against KEGG ---
set.seed(42)
gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  minGSSize    = 5,
  maxGSSize    = 500,
  pvalueCutoff = 1,     # keep everything, we'll filter after
  eps          = 0,
  verbose      = FALSE
)

gsea_df <- as.data.frame(gsea_kegg)
cat("\nTotal KEGG pathways tested via GSEA:", nrow(gsea_df), "\n")

# --- Filter to IVD-relevant pathways ---
gsea_ivd_df <- gsea_df[gsea_df$ID %in% ivd_relevant_pathways, ]
gsea_ivd_df <- gsea_ivd_df[order(gsea_ivd_df$pvalue), ]

cat("\nIVD-relevant pathways found via GSEA:\n")
print(gsea_ivd_df[, c("ID", "Description", "NES", "pvalue", "p.adjust", "setSize")])

write.csv(gsea_ivd_df, "GSEA_KEGG_IVD_relevant_sNPC_vs_sNPCMG.csv", row.names = FALSE)

# ============================================================
# Top 10 IVD pathways - bar graph with NES direction
# ============================================================
library(ggplot2)
library(dplyr)

top10_gsea_ivd <- gsea_ivd_df %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

top10_gsea_ivd$Direction <- ifelse(top10_gsea_ivd$NES > 0, "Up in sNPC-MG", "Up in sNPC")

top10_gsea_ivd$Description <- factor(
  top10_gsea_ivd$Description,
  levels = top10_gsea_ivd$Description[order(top10_gsea_ivd$NES)]
)

p_gsea_ivd <- ggplot(top10_gsea_ivd, aes(x = Description, y = NES, fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Up in sNPC" = "#4575b4", "Up in sNPC-MG" = "#d73027")) +
  labs(title = "Top 10 IVD-Relevant KEGG Pathways (GSEA): sNPC vs sNPC-MG",
       x = "", y = "Normalized Enrichment Score (NES)", fill = "Direction") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

print(p_gsea_ivd)

ggsave("Top10_IVD_KEGG_GSEA_barplot.pdf", plot = p_gsea_ivd, width = 9, height = 6)
cat("\nSaved: Top10_IVD_KEGG_GSEA_barplot.pdf\n")


############################### senecent marker ###############


# ============================================================
# Senescence Marker Panel Analysis: sNPC vs sNPC-MG
# ============================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# --- Curated senescence marker panel ---
# Cell cycle arrest / senescence effectors, SASP factors, structural markers
senescence_markers <- c(
  "CDKN1A",   # p21
  "CDKN2A",   # p16INK4a
  "TP53",     # p53
  "GLB1",     # senescence-associated beta-galactosidase (lysosomal)
  "SERPINE1", # PAI-1, classic SASP factor
  "IL6",      # SASP
  "IL1B",     # SASP
  "IL8",      # SASP / CXCL8
  "CXCL8",    # alt symbol
  "MMP3",     # SASP-associated matrix remodeling
  "MMP12",    # SASP
  "LMNB1",    # Lamin B1 -- DOWNregulated in senescence
  "HMGB1",    # released in senescence
  "IGFBP3",   # SASP
  "IGFBP7",   # SASP
  "CCL2",     # SASP (MCP-1)
  "TNF",      # SASP
  "RB1",      # retinoblastoma, senescence pathway
  "TERT"      # telomerase (often down in senescence)
)

# --- Match against log_mat_imp rownames (case-insensitive) ---
rn <- rownames(log_mat_imp)
rn_upper <- toupper(trimws(rn))
markers_upper <- toupper(senescence_markers)

matched_idx <- match(markers_upper, rn_upper)
found <- !is.na(matched_idx)

matched_rows <- rn[matched_idx[found]]
matched_query <- senescence_markers[found]

cat("Senescence markers found in data:\n")
print(data.frame(query = matched_query, matched_rowname = matched_rows))

missing_markers <- senescence_markers[!found]
if (length(missing_markers) > 0) {
  cat("\nNOT found in log_mat_imp:", paste(missing_markers, collapse = ", "), "\n")
  cat("(Will check larger unfiltered matrices below)\n")
}

# --- Check unfiltered matrix for any markers missing from log_mat_imp ---
if (length(missing_markers) > 0 && exists("log_mat")) {
  rn2 <- toupper(trimws(rownames(log_mat)))
  found_in_full <- missing_markers[toupper(missing_markers) %in% rn2]
  if (length(found_in_full) > 0) {
    cat("Found in unfiltered log_mat (but filtered out of log_mat_imp):", paste(found_in_full, collapse = ", "), "\n")
  }
}

matched_unique <- unique(matched_rows)
if (length(matched_unique) == 0) stop("No senescence markers found. Check gene symbol formatting.")

# ============================================================
# Subset to sNPC and sNPC-MG samples only
# ============================================================
keep_samples_sen <- sample_metadata$Sample[sample_metadata$Group %in% c("sNPC", "sNPC-MG")]
keep_samples_sen <- intersect(keep_samples_sen, colnames(log_mat_imp))

sub_mat_sen <- log_mat_imp[matched_unique, keep_samples_sen, drop = FALSE]

# --- Long format ---
data_long_sen <- melt(as.matrix(sub_mat_sen), varnames = c("Marker", "Sample"), value.name = "log2_Intensity")
data_long_sen <- merge(data_long_sen, sample_metadata[, c("Sample", "Group")], by = "Sample")
data_long_sen <- data_long_sen[!is.na(data_long_sen$log2_Intensity), ]
data_long_sen$Group <- factor(data_long_sen$Group, levels = c("sNPC", "sNPC-MG"))

# ============================================================
# Stats: mean, SD, t-test per marker
# ============================================================
summary_sen <- data_long_sen %>%
  group_by(Marker, Group) %>%
  summarise(Mean = mean(log2_Intensity), SD = sd(log2_Intensity), N = n(), .groups = "drop")

ttest_results <- data_long_sen %>%
  group_by(Marker) %>%
  summarise(
    p_value = tryCatch(t.test(log2_Intensity ~ Group)$p.value, error = function(e) NA),
    .groups = "drop"
  )

# Pull matching logFC/adj.P from your existing limma results if available
if (exists("de_results")) {
  de_sub <- de_results[rownames(de_results) %in% matched_unique, c("logFC", "P.Value", "adj.P.Val")]
  de_sub$Marker <- rownames(de_sub)
  ttest_results <- merge(ttest_results, de_sub, by = "Marker", all.x = TRUE)
}

ttest_results <- ttest_results[order(ttest_results$p_value), ]
cat("\nSenescence marker statistics (sNPC vs sNPC-MG):\n")
print(ttest_results)

write.csv(ttest_results, "Senescence_markers_stats_sNPC_vs_sNPCMG.csv", row.names = FALSE)

# ============================================================
# Bar graph: Mean +/- SD per marker, colored by group
# ============================================================
group_colors_sen <- c("sNPC" = "#4575b4", "sNPC-MG" = "#d73027")

p_sen_bar <- ggplot(summary_sen, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.25, linewidth = 0.7) +
  facet_wrap(~ Marker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors_sen) +
  labs(title = "Senescence Markers: sNPC vs sNPC-MG (Mean ± SD)",
       x = "", y = "log2 Intensity") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

print(p_sen_bar)
ggsave("Senescence_markers_barplot.pdf", plot = p_sen_bar, width = 11, height = 8)

# ============================================================
# Violin plot: distribution per marker
# ============================================================
p_sen_violin <- ggplot(data_long_sen, aes(x = Group, y = log2_Intensity, fill = Group)) +
  geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.6) +
  geom_point(data = summary_sen, aes(x = Group, y = Mean),
             shape = 23, size = 2.5, fill = "white", color = "black", inherit.aes = FALSE) +
  facet_wrap(~ Marker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors_sen) +
  labs(title = "Senescence Marker Distribution: sNPC vs sNPC-MG",
       x = "", y = "log2 Intensity") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

print(p_sen_violin)
ggsave("Senescence_markers_violinplot.pdf", plot = p_sen_violin, width = 11, height = 8)

cat("\nSaved: Senescence_markers_barplot.pdf, Senescence_markers_violinplot.pdf, Senescence_markers_stats_sNPC_vs_sNPCMG.csv\n")

#################################### violine  sNPC vs sNPC-MG###########################

# ============================================================
# General Violin Plot: sNPC vs sNPC-MG (Top DE proteins, whole proteome)
# ============================================================
library(ggplot2)
library(reshape2)
library(dplyr)

# --- Use sig_genes (from limma DE results, whole proteome) ---
cat("Total significant DE proteins (adj.P < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")

# --- Select top proteins to visualize (by adj.P.Val, capped at ~20 for readability) ---
top_n <- 20
top_proteins <- sig_genes %>%
  arrange(adj.P.Val) %>%
  slice_head(n = top_n) %>%
  pull(Protein)

cat("\nTop", length(top_proteins), "DE proteins selected for plotting:\n")
print(top_proteins)

# ============================================================
# Subset expression matrix to top proteins, sNPC + sNPC-MG samples
# ============================================================
keep_samples_all <- sample_metadata$$Sample[sample_metadata$$Group %in% c("sNPC", "sNPC-MG")]
keep_samples_all <- intersect(keep_samples_all, colnames(log_mat_imp))

sub_mat_all <- log_mat_imp[top_proteins, keep_samples_all, drop = FALSE]

data_long_all <- melt(as.matrix(sub_mat_all), varnames = c("Protein", "Sample"), value.name = "log2_Intensity")
data_long_all <- merge(data_long_all, sample_metadata[, c("Sample", "Group")], by = "Sample")
data_long_all <- data_long_all[!is.na(data_long_all$log2_Intensity), ]
data_long_all$$Group <- factor(data_long_all$$Group, levels = c("sNPC", "sNPC-MG"))

# ============================================================
# Direction labeling (Up/Down in sNPC-MG)
# ============================================================
sig_genes$$Direction <- ifelse(sig_genes$$logFC > 0, "Upregulated (sNPC-MG)", "Downregulated (sNPC-MG)")

up_proteins   <- sig_genes$$Protein[sig_genes$$Direction == "Upregulated (sNPC-MG)"]
down_proteins <- sig_genes$$Protein[sig_genes$$Direction == "Downregulated (sNPC-MG)"]

cat("\n=== Upregulated in sNPC-MG ===\n")
print(up_proteins)

cat("\n=== Downregulated in sNPC-MG ===\n")
print(down_proteins)

# --- Build facet labels with direction arrows ---
dir_lookup_all <- setNames(sig_genes$$Direction, sig_genes$$Protein)

label_map_all <- sapply(top_proteins, function(p) {
  d <- dir_lookup_all[p]
  ifelse(d == "Upregulated (sNPC-MG)", paste0(p, " \u2191"), paste0(p, " \u2193"))
})

data_long_all$$Protein_labeled <- factor(data_long_all$$Protein,
                                         levels = names(label_map_all),
                                         labels = label_map_all)

# ============================================================
# Violin plot
# ============================================================
group_colors <- c("sNPC" = "#4575b4", "sNPC-MG" = "#d73027")

p_violin_all <- ggplot(data_long_all, aes(x = Group, y = log2_Intensity, fill = Group)) +
  geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
  geom_jitter(width = 0.1
              
              
              
              keep_samples_all <- sample_metadata$Sample[sample_metadata$Group %in% c("sNPC", "sNPC-MG")]
              
              
              
              sub_mat_all <- log_mat_imp[top_proteins, keep_samples_all, drop = FALSE]
              
              data_long_all <- melt(as.matrix(sub_mat_all), varnames = c("Protein", "Sample"), value.name = "log2_Intensity")
              data_long_all <- merge(data_long_all, sample_metadata[, c("Sample", "Group")], by = "Sample")
              data_long_all <- data_long_all[!is.na(data_long_all$log2_Intensity), ]
              data_long_all$Group <- factor(data_long_all$Group, levels = c("sNPC", "sNPC-MG"))
              
              print(keep_samples_all)
              
              
              sub_mat_all <- log_mat_imp[top_proteins, keep_samples_all, drop = FALSE]
              
              data_long_all <- melt(as.matrix(sub_mat_all), varnames = c("Protein", "Sample"), value.name = "log2_Intensity")
              data_long_all <- merge(data_long_all, sample_metadata[, c("Sample", "Group")], by = "Sample")
              data_long_all <- data_long_all[!is.na(data_long_all$log2_Intensity), ]
              data_long_all$Group <- factor(data_long_all$Group, levels = c("sNPC", "sNPC-MG"))
              
              # --- Direction labeling ---
              sig_genes$Direction <- ifelse(sig_genes$logFC > 0, "Upregulated (sNPC-MG)", "Downregulated (sNPC-MG)")
              
              up_proteins   <- sig_genes$Protein[sig_genes$Direction == "Upregulated (sNPC-MG)"]
              down_proteins <- sig_genes$Protein[sig_genes$Direction == "Downregulated (sNPC-MG)"]
              
              cat("\n=== Upregulated in sNPC-MG ===\n")
              print(up_proteins)
              
              cat("\n=== Downregulated in sNPC-MG ===\n")
              print(down_proteins)
              
              # --- Facet labels ---
              dir_lookup_all <- setNames(sig_genes$Direction, sig_genes$Protein)
              
              label_map_all <- sapply(top_proteins, function(p) {
                d <- dir_lookup_all[p]
                ifelse(d == "Upregulated (sNPC-MG)", paste0(p, " \u2191"), paste0(p, " \u2193"))
              })
              
              data_long_all$Protein_labeled <- factor(data_long_all$Protein,
                                                      levels = names(label_map_all),
                                                      labels = label_map_all)
              
              # --- Violin plot ---
              group_colors <- c("sNPC" = "#4575b4", "sNPC-MG" = "#d73027")
              
              p_violin_all <- ggplot(data_long_all, aes(x = Group, y = log2_Intensity, fill = Group)) +
                geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
                geom_jitter(width = 0.1, size = 1.2, alpha = 0.6, color = "black") +
                stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", color = "black") +
                facet_wrap(~ Protein_labeled, scales = "free_y", ncol = 5) +
                scale_fill_manual(values = group_colors) +
                labs(title = "Top Differentially Expressed Proteins: sNPC vs sNPC-MG",
                     subtitle = "\u2191 = Upregulated in sNPC-MG | \u2193 = Downregulated in sNPC-MG",
                     x = "", y = "log2 Intensity") +
                theme_bw(base_size = 10.5) +
                theme(
                  axis.text.x = element_text(face = "bold"),
                  strip.text = element_text(face = "bold", size = 9),
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
                  legend.position = "none"
                )
              
              print(p_violin_all)
              
              ggsave("Top_DE_proteins_violinplot_sNPC_vs_sNPCMG.pdf", plot = p_violin_all, width = 13, height = 9)
              cat("\nViolin plot saved: Top_DE_proteins_violinplot_sNPC_vs_sNPCMG.pdf\n")
              ################### VOCONO PLOT#############
              
              
              # ============================================================
              # Volcano Plot: sNPC vs sNPC-MG
              # ============================================================
              library(ggplot2)
              library(ggrepel)
              library(dplyr)
              
              # --- Thresholds ---
              fc_cutoff <- 1        # |log2FC| > 1
              p_cutoff  <- 0.05      # adj.P.Val < 0.05
              
              # --- Prepare data ---
              volcano_df <- de_results
              volcano_df$Protein <- rownames(volcano_df)
              
              volcano_df$Category <- case_when(
                volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              volcano_df$Category <- factor(volcano_df$Category,
                                            levels = c("Upregulated (sNPC-MG)", "Downregulated (sNPC-MG)", "Not significant"))
              
              cat("Category counts:\n")
              print(table(volcano_df$Category))
              
              # --- Select top proteins to label (by adj.P.Val, within significant categories) ---
              top_label_n <- 10
              
              top_up <- volcano_df %>%
                filter(Category == "Upregulated (sNPC-MG)") %>%
                arrange(adj.P.Val) %>%
                slice_head(n = top_label_n)
              
              top_down <- volcano_df %>%
                filter(Category == "Downregulated (sNPC-MG)") %>%
                arrange(adj.P.Val) %>%
                slice_head(n = top_label_n)
              
              labels_df <- bind_rows(top_up, top_down)
              
              cat("\nTop upregulated proteins (sNPC-MG):\n")
              print(top_up$Protein)
              
              cat("\nTop downregulated proteins (sNPC-MG):\n")
              print(top_down$Protein)
              
              # ============================================================
              # Volcano plot
              # ============================================================
              volcano_colors <- c(
                "Upregulated (sNPC-MG)"   = "#d73027",
                "Downregulated (sNPC-MG)" = "#4575b4",
                "Not significant"         = "grey70"
              )
              
              p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
                geom_point(alpha = 0.7, size = 1.8) +
                scale_color_manual(values = volcano_colors) +
                geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
                geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey40") +
                geom_text_repel(data = labels_df, aes(label = Protein),
                                size = 3.2, max.overlaps = 20, box.padding = 0.4,
                                segment.color = "grey50", color = "black", show.legend = FALSE) +
                labs(title = "Volcano Plot: sNPC vs sNPC-MG",
                     subtitle = paste0("Cutoffs: |log2FC| > ", fc_cutoff, ", adj.P < ", p_cutoff),
                     x = "log2 Fold Change (sNPC-MG vs sNPC)",
                     y = "-log10(adjusted p-value)",
                     color = "") +
                theme_bw(base_size = 12) +
                theme(
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
                  legend.position = "bottom"
                )
              
              print(p_volcano)
              
              ggsave("Volcano_sNPC_vs_sNPCMG.pdf", plot = p_volcano, width = 8, height = 7)
              cat("\nVolcano plot saved: Volcano_sNPC_vs_sNPCMG.pdf\n")
              
              write.csv(volcano_df[, c("Protein", "logFC", "P.Value", "adj.P.Val", "Category")],
                        "Volcano_data_sNPC_vs_sNPCMG.csv", row.names = FALSE)
              
              
              
              # ============================================================
              # Volcano Plot: sNPC vs sNPC-MG — Publication Style
              # ============================================================
              library(ggplot2)
              library(ggrepel)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              volcano_df <- de_results
              volcano_df$Protein <- rownames(volcano_df)
              
              volcano_df$Category <- case_when(
                volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              volcano_df$Category <- factor(volcano_df$Category,
                                            levels = c("Upregulated (sNPC-MG)", "Downregulated (sNPC-MG)", "Not significant"))
              
              top_label_n <- 10
              
              top_up <- volcano_df %>%
                filter(Category == "Upregulated (sNPC-MG)") %>%
                arrange(adj.P.Val) %>%
                slice_head(n = top_label_n)
              
              top_down <- volcano_df %>%
                filter(Category == "Downregulated (sNPC-MG)") %>%
                arrange(adj.P.Val) %>%
                slice_head(n = top_label_n)
              
              labels_df <- bind_rows(top_up, top_down)
              
              # ============================================================
              # Publication-style volcano plot
              # ============================================================
              volcano_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4",
                "Not significant"         = "grey75"
              )
              
              p_volcano_pub <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
                geom_point(alpha = 0.8, size = 2) +
                scale_color_manual(values = volcano_colors) +
                geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
                geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
                geom_text_repel(
                  data = labels_df, aes(label = Protein),
                  size = 3.6, fontface = "bold", color = "black",
                  max.overlaps = 30, box.padding = 0.45, point.padding = 0.3,
                  segment.color = "grey30", segment.size = 0.4,
                  show.legend = FALSE
                ) +
                labs(
                  title = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = expression(bold(paste("-log"[10], "(adjusted ", italic("p"), "-value)"))),
                  color = ""
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 13),
                  axis.text = element_text(face = "bold", size = 11, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(face = "bold", size = 11),
                  legend.title = element_blank(),
                  legend.position = "top",
                  legend.key = element_blank(),
                  panel.grid = element_blank(),
                  plot.margin = margin(15, 15, 10, 10)
                ) +
                guides(color = guide_legend(override.aes = list(size = 3)))
              
              print(p_volcano_pub)
              
              ggsave("Volcano_sNPC_vs_sNPCMG_publication.pdf", plot = p_volcano_pub, width = 7.5, height = 7, dpi = 600)
              ggsave("Volcano_sNPC_vs_sNPCMG_publication.png", plot = p_volcano_pub, width = 7.5, height = 7, dpi = 600)
              
              cat("Saved: Volcano_sNPC_vs_sNPCMG_publication.pdf / .png (600 dpi)\n")
      
              
              # --- Upregulated labels: repel to the right ---
              geom_text_repel(
                data = top_up, aes(label = Protein),
                size = 3.6, fontface = "bold", color = "black",
                nudge_x = top_up$nudge_x_val,
                nudge_y = top_up$nudge_y_val,
                direction = "both",
                hjust = 0,
                segment.color = "grey30", segment.size = 0.4,
                min.segment.length = 0,
                box.padding = 0.8,
                point.padding = 0.4,
                force = 15,
                force_pull = 0.1,
                max.iter = 50000,
                max.time = 3,
                max.overlaps = Inf,
                seed = 42,
                show.legend = FALSE
              ) +
                
                # --- Downregulated labels: repel to the left ---
                geom_text_repel(
                  data = top_down, aes(label = Protein),
                  size = 3.6, fontface = "bold", color = "black",
                  nudge_x = top_down$nudge_x_val,
                  nudge_y = top_down$nudge_y_val,
                  direction = "both",
                  hjust = 1,
                  segment.color = "grey30", segment.size = 0.4,
                  min.segment.length = 0,
                  box.padding = 0.8,
                  point.padding = 0.4,
                  force = 15,
                  force_pull = 0.1,
                  max.iter = 50000,
                  max.time = 3,
                  max.overlaps = Inf,
                  seed = 42,
                  show.legend = FALSE
                ) +
                
 ############################ ranked bar plot ##########################################
              
              
              # ============================================================
              # Ranked Bar Plot: sNPC vs sNPC-MG
              # Top Up/Down proteins ranked by log2FC
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              bar_df <- de_results
              bar_df$Protein <- rownames(bar_df)
              
              bar_df$Category <- case_when(
                bar_df$adj.P.Val < p_cutoff & bar_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                bar_df$adj.P.Val < p_cutoff & bar_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              # --- Select top N up and top N down proteins ---
              top_n <- 10
              
              top_up <- bar_df %>%
                filter(Category == "Upregulated (sNPC-MG)") %>%
                arrange(desc(logFC)) %>%
                slice_head(n = top_n)
              
              top_down <- bar_df %>%
                filter(Category == "Downregulated (sNPC-MG)") %>%
                arrange(logFC) %>%
                slice_head(n = top_n)
              
              plot_df <- bind_rows(top_up, top_down)
              
              # --- Order proteins by logFC for ranked display ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Top upregulated proteins (sNPC-MG):\n")
              print(top_up$Protein)
              
              cat("\nTop downregulated proteins (sNPC-MG):\n")
              print(top_down$Protein)
              
              # ============================================================
              # Ranked bar plot
              # ============================================================
              bar_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_bar <- ggplot(plot_df, aes(x = Protein, y = logFC, fill = Category)) +
                geom_col(width = 0.7, color = "black", linewidth = 0.3) +
                coord_flip() +
                scale_fill_manual(values = bar_colors) +
                labs(
                  title = "Top Differentially Expressed Proteins",
                  subtitle = "sNPC vs sNPC-MG",
                  x = "",
                  y = expression(bold(paste("log"[2], " Fold Change"))),
                  fill = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(face = "bold", size = 11),
                  legend.title = element_blank(),
                  legend.position = "top",
                  panel.grid = element_blank()
                )
              
              print(p_bar)
              
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.pdf", plot = p_bar, width = 7, height = 8, dpi = 600)
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.png", plot = p_bar, width = 7, height = 8, dpi = 600)
              
              cat("\nSaved: RankedBarPlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              
              class(p_bar)
              print(p_bar)
              
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.pdf", plot = p_bar, width = 7, height = 8, dpi = 600)
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.png", plot = p_bar, width = 7, height = 8, dpi = 600)
              
              
              # ============================================================
              # Ranked Bar Plot: sNPC vs sNPC-MG
              # Top Up/Down proteins ranked by log2FC
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              bar_df <- de_results
              bar_df$Protein <- rownames(bar_df)
              
              bar_df$Category <- case_when(
                bar_df$adj.P.Val < p_cutoff & bar_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                bar_df$adj.P.Val < p_cutoff & bar_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              # --- Select top N up and top N down proteins ---
              top_n <- 10
              
              top_up <- bar_df %>%
                filter(Category == "Upregulated (sNPC-MG)") %>%
                arrange(desc(logFC)) %>%
                slice_head(n = top_n)
              
              top_down <- bar_df %>%
                filter(Category == "Downregulated (sNPC-MG)") %>%
                arrange(logFC) %>%
                slice_head(n = top_n)
              
              plot_df <- bind_rows(top_up, top_down)
              
              # --- Order proteins by logFC for ranked display ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Top upregulated proteins (sNPC-MG):\n")
              print(top_up$Protein)
              
              cat("\nTop downregulated proteins (sNPC-MG):\n")
              print(top_down$Protein)
              
              # ============================================================
              # Ranked bar plot
              # ============================================================
              bar_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_bar <- ggplot(plot_df, aes(x = Protein, y = logFC, fill = Category)) +
                geom_col(width = 0.7, color = "black", linewidth = 0.3) +
                coord_flip() +
                scale_fill_manual(values = bar_colors) +
                labs(
                  title = "Top Differentially Expressed Proteins",
                  subtitle = "sNPC vs sNPC-MG",
                  x = "",
                  y = expression(bold(paste("log"[2], " Fold Change"))),
                  fill = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(face = "bold", size = 11),
                  legend.title = element_blank(),
                  legend.position = "top",
                  panel.grid = element_blank()
                )
              
              print(p_bar)
              
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.pdf", plot = p_bar, width = 7, height = 8, dpi = 600)
              ggsave("RankedBarPlot_sNPC_vs_sNPCMG.png", plot = p_bar, width = 7, height = 8, dpi = 600)
              
              cat("\nSaved: RankedBarPlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              ################### ranked bar plot mitondrial proitein############
              
              
              # ============================================================
              # Ranked Bar Plot: Mitochondrial Proteins Only
              # sNPC vs sNPC-MG
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              mito_df <- de_results
              mito_df$Protein <- rownames(mito_df)
              
              mito_df$Category <- case_when(
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              # ============================================================
              # OPTION 1: You have a vector of mitochondrial protein/gene names
              # Edit this list to match your mitochondrial proteins of interest
              # ============================================================
              mito_genes <- c("NDUFA1", "NDUFS1", "SDHA", "SDHB", "UQCRC1", "UQCRC2",
                              "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B", "TFAM", "MFN1",
                              "MFN2", "OPA1", "DNM1L", "PINK1", "PRKN", "TOMM20", "VDAC1")
              
              mito_df_filtered <- mito_df %>%
                filter(Protein %in% mito_genes)
              
              # ============================================================
              # OPTION 2 (alternative): If your de_results already has a column
              # flagging mitochondrial proteins (e.g., "Is_Mitochondrial" = TRUE/FALSE),
              # comment out Option 1 above and uncomment this instead:
              # ============================================================
              # mito_df_filtered <- mito_df %>% filter(Is_Mitochondrial == TRUE)
              
              # --- Keep only significant mitochondrial proteins for the plot ---
              plot_df <- mito_df_filtered %>%
                filter(Category != "Not significant") %>%
                arrange(logFC)
              
              # --- Order proteins by logFC for ranked display ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Mitochondrial proteins included in plot:\n")
              print(plot_df$Protein)
              
              # ============================================================
              # Ranked bar plot
              # ============================================================
              bar_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_mito_bar <- ggplot(plot_df, aes(x = Protein, y = logFC, fill = Category)) +
                geom_col(width = 0.7, color = "black", linewidth = 0.3) +
                coord_flip() +
                scale_fill_manual(values = bar_colors) +
                labs(
                  title = "Mitochondrial Protein Expression",
                  subtitle = "sNPC vs sNPC-MG",
                  x = "",
                  y = expression(bold(paste("log"[2], " Fold Change"))),
                  fill = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(face = "bold", size = 11),
                  legend.title = element_blank(),
                  legend.position = "top",
                  panel.grid = element_blank()
                )
              
              print(p_mito_bar)
              
              ggsave("MitoRankedBarPlot_sNPC_vs_sNPCMG.pdf", plot = p_mito_bar, width = 7, height = 8, dpi = 600)
              ggsave("MitoRankedBarPlot_sNPC_vs_sNPCMG.png", plot = p_mito_bar, width = 7, height = 8, dpi = 600)
              
              cat("\nSaved: MitoRankedBarPlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              ############################ dot plot##################################
              
              
              # ============================================================
              # Bubble Plot: Mitochondrial Proteins
              # sNPC vs sNPC-MG
              # X-axis = log2FC | Dot size = -log10(p-value) | Color = direction
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              mito_df <- de_results
              mito_df$Protein <- rownames(mito_df)
              
              mito_df$Category <- case_when(
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              mito_df$negLogP <- -log10(mito_df$adj.P.Val)
              
              # ============================================================
              # OPTION 1: You have a vector of mitochondrial protein/gene names
              # Edit this list to match your mitochondrial proteins of interest
              # ============================================================
              mito_genes <- c("NDUFA1", "NDUFS1", "SDHA", "SDHB", "UQCRC1", "UQCRC2",
                              "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B", "TFAM", "MFN1",
                              "MFN2", "OPA1", "DNM1L", "PINK1", "PRKN", "TOMM20", "VDAC1")
              
              mito_df_filtered <- mito_df %>%
                filter(Protein %in% mito_genes)
              
              # ============================================================
              # OPTION 2 (alternative): If de_results already has a column flagging
              # mitochondrial proteins (e.g., "Is_Mitochondrial" = TRUE/FALSE),
              # comment out Option 1 above and uncomment this instead:
              # ============================================================
              # mito_df_filtered <- mito_df %>% filter(Is_Mitochondrial == TRUE)
              
              # --- Keep only significant mitochondrial proteins ---
              plot_df <- mito_df_filtered %>%
                filter(Category != "Not significant") %>%
                arrange(logFC)
              
              # --- Order proteins by logFC for ranked display on y-axis ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Mitochondrial proteins included in plot:\n")
              print(plot_df$Protein)
              
              # ============================================================
              # Bubble plot
              # ============================================================
              bubble_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_bubble <- ggplot(plot_df, aes(x = logFC, y = Protein, size = negLogP, color = Category)) +
                geom_point(alpha = 0.85) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
                scale_color_manual(values = bubble_colors) +
                scale_size_continuous(name = expression(bold(paste("-log"[10], "(adj. p-value)"))),
                                      range = c(3, 10)) +
                labs(
                  title = "Mitochondrial Protein Expression",
                  subtitle = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = "",
                  color = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(face = "bold", size = 10),
                  legend.position = "right",
                  panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                guides(color = guide_legend(override.aes = list(size = 5)))
              
              print(p_bubble)
              
              ggsave("MitoBubblePlot_sNPC_vs_sNPCMG.pdf", plot = p_bubble, width = 8, height = 8, dpi = 600)
              ggsave("MitoBubblePlot_sNPC_vs_sNPCMG.png", plot = p_bubble, width = 8, height = 8, dpi = 600)
              
              cat("\nSaved: MitoBubblePlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              colnames(de_results)
              # ============================================================
              # Bubble Plot: Mitochondrial Proteins
              # sNPC vs sNPC-MG
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              mito_df <- de_results
              mito_df$Protein <- rownames(mito_df)
              
              mito_df$Category <- case_when(
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                mito_df$adj.P.Val < p_cutoff & mito_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              mito_genes <- c("NDUFA1", "NDUFS1", "SDHA", "SDHB", "UQCRC1", "UQCRC2",
                              "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B", "TFAM", "MFN1",
                              "MFN2", "OPA1", "DNM1L", "PINK1", "PRKN", "TOMM20", "VDAC1")
              
              mito_df_filtered <- mito_df %>%
                filter(Protein %in% mito_genes)
              
              plot_df <- mito_df_filtered %>%
                filter(Category != "Not significant") %>%
                mutate(negLogP = -log10(adj.P.Val)) %>%   # recompute here to guarantee it exists
                arrange(logFC)
              
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Columns in plot_df:\n")
              print(colnames(plot_df))
              
              cat("\nMitochondrial proteins included in plot:\n")
              print(plot_df$Protein)
              
              # ============================================================
              # Bubble plot
              # ============================================================
              bubble_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_bubble <- ggplot(plot_df, aes(x = logFC, y = Protein, size = negLogP, color = Category)) +
                geom_point(alpha = 0.85) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
                scale_color_manual(values = bubble_colors) +
                scale_size_continuous(name = expression(bold(paste("-log"[10], "(adj. p-value)"))),
                                      range = c(3, 10)) +
                labs(
                  title = "Mitochondrial Protein Expression",
                  subtitle = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = "",
                  color = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(face = "bold", size = 10),
                  legend.position = "right",
                  panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                guides(color = guide_legend(override.aes = list(size = 5)))
              
              print(p_bubble)
              
              ggsave("MitoBubblePlot_sNPC_vs_sNPCMG.pdf", plot = p_bubble, width = 8, height = 8, dpi = 600)
              ggsave("MitoBubblePlot_sNPC_vs_sNPCMG.png", plot = p_bubble, width = 8, height = 8, dpi = 600)
              
              cat("\nSaved: MitoBubblePlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              
              ##################### PROINFLMMATORY##########################
              
              
              # ============================================================
              # Bubble Plot: Top 10 Proinflammatory Proteins
              # sNPC vs sNPC-MG
              # X-axis = log2FC | Dot size = -log10(p-value) | Color = direction
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              inflam_df <- de_results
              inflam_df$Protein <- rownames(inflam_df)
              
              inflam_df$Category <- case_when(
                inflam_df$adj.P.Val < p_cutoff & inflam_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                inflam_df$adj.P.Val < p_cutoff & inflam_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              # ============================================================
              # Proinflammatory gene/protein list — EDIT to match your rownames
              # ============================================================
              proinflam_genes <- c("IL1B", "IL6", "TNF", "NFKB1", "NLRP3", "CXCL8",
                                   "CCL2", "IL18", "PTGS2", "IL1A", "CXCL10", "CCL5",
                                   "IL12A", "IFNG", "TLR4", "MYD88", "STAT1", "NOS2",
                                   "CASP1", "IL23A")
              
              inflam_df_filtered <- inflam_df %>%
                filter(Protein %in% proinflam_genes)
              
              # --- Select top 10 proinflammatory proteins by significance ---
              top_n <- 10
              
              plot_df <- inflam_df_filtered %>%
                filter(Category != "Not significant") %>%
                mutate(negLogP = -log10(adj.P.Val)) %>%
                arrange(adj.P.Val) %>%      # rank by significance
                slice_head(n = top_n) %>%
                arrange(logFC)
              
              # --- Order proteins by logFC for y-axis display ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Top 10 proinflammatory proteins included in plot:\n")
              print(plot_df$Protein)
              
              # ============================================================
              # Bubble plot
              # ============================================================
              bubble_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              p_bubble_inflam <- ggplot(plot_df, aes(x = logFC, y = Protein, size = negLogP, color = Category)) +
                geom_point(alpha = 0.85) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
                scale_color_manual(values = bubble_colors) +
                scale_size_continuous(name = expression(bold(paste("-log"[10], "(adj. p-value)"))),
                                      range = c(3, 10)) +
                labs(
                  title = "Top 10 Proinflammatory Protein Enrichment",
                  subtitle = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = "",
                  color = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(face = "bold", size = 10),
                  legend.position = "right",
                  panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                guides(color = guide_legend(override.aes = list(size = 5)))
              
              print(p_bubble_inflam)
              
              ggsave("ProinflammatoryBubblePlot_sNPC_vs_sNPCMG.pdf", plot = p_bubble_inflam, width = 8, height = 7, dpi = 600)
              ggsave("ProinflammatoryBubblePlot_sNPC_vs_sNPCMG.png", plot = p_bubble_inflam, width = 8, height = 7, dpi = 600)
              
              cat("\nSaved: ProinflammatoryBubblePlot_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              
              # How many of your proinflammatory genes actually exist in the data at all?
              sum(proinflam_genes %in% rownames(de_results))
              
              # What do the actual row names look like?
              head(rownames(de_results), 20)
              
              
              
              # How many of your listed proinflammatory genes are actually present?
              matched <- proinflam_genes[proinflam_genes %in% rownames(de_results)]
              matched
              
              # How many matched?
              length(matched)
              
              
              proinflam_genes_expanded <- c(
                # NF-kB / signaling
                "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "CHUK", "MYD88", "TRAF6", "TRAF2",
                # JAK-STAT
                "STAT1", "STAT3", "JAK1", "JAK2", "SOCS3",
                # Cytokines/chemokines
                "IL1B", "IL6", "TNF", "CXCL8", "CCL2", "CCL5", "CXCL10", "IL18", "IL1A", "IL33",
                # Prostaglandin / COX pathway
                "PTGS2", "PTGES", "ALOX5",
                # Inflammasome
                "NLRP3", "CASP1", "PYCARD", "IL1RN",
                # TLR pathway
                "TLR2", "TLR4", "TLR3", "TICAM1",
                # Acute phase / complement
                "CRP", "SAA1", "SAA2", "C3", "C1QA", "C1QB", "C1QC", "CFB",
                # S100 alarmins
                "S100A8", "S100A9", "S100B",
                # Other inflammation-associated
                "TNFAIP3", "TNFRSF1A", "TNFRSF1B", "ICAM1", "VCAM1", "PTPN22", "CD68", "AIF1"
              )
              
              # Check overlap with your dataset
              matched_expanded <- proinflam_genes_expanded[proinflam_genes_expanded %in% rownames(de_results)]
              matched_expanded
              length(matched_expanded)
              
              
              
              # ============================================================
              # Bubble Plot: All Significant Proinflammatory Proteins
              # sNPC vs sNPC-MG
              # X-axis = log2FC | Dot size = -log10(p-value) | Color = direction
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              fc_cutoff <- 1
              p_cutoff  <- 0.05
              
              # --- Prepare DE data ---
              inflam_df <- de_results
              inflam_df$Protein <- rownames(inflam_df)
              
              inflam_df$Category <- case_when(
                inflam_df$adj.P.Val < p_cutoff & inflam_df$logFC >  fc_cutoff ~ "Upregulated (sNPC-MG)",
                inflam_df$adj.P.Val < p_cutoff & inflam_df$logFC < -fc_cutoff ~ "Downregulated (sNPC-MG)",
                TRUE ~ "Not significant"
              )
              
              # ============================================================
              # Expanded proinflammatory gene panel
              # ============================================================
              proinflam_genes_expanded <- c(
                "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "CHUK", "MYD88", "TRAF6", "TRAF2",
                "STAT1", "STAT3", "JAK1", "JAK2", "SOCS3",
                "IL1B", "IL6", "TNF", "CXCL8", "CCL2", "CCL5", "CXCL10", "IL18", "IL1A", "IL33",
                "PTGS2", "PTGES", "ALOX5",
                "NLRP3", "CASP1", "PYCARD", "IL1RN",
                "TLR2", "TLR4", "TLR3", "TICAM1",
                "CRP", "SAA1", "SAA2", "C3", "C1QA", "C1QB", "C1QC", "CFB",
                "S100A8", "S100A9", "S100B",
                "TNFAIP3", "TNFRSF1A", "TNFRSF1B", "ICAM1", "VCAM1", "PTPN22", "CD68", "AIF1"
              )
              
              inflam_df_filtered <- inflam_df %>%
                filter(Protein %in% proinflam_genes_expanded)
              
              # --- Keep ALL significant proinflammatory proteins (no top-N cutoff) ---
              plot_df <- inflam_df_filtered %>%
                filter(Category != "Not significant") %>%
                mutate(negLogP = -log10(adj.P.Val)) %>%
                arrange(logFC)
              
              # --- Order proteins by logFC for y-axis display ---
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              cat("Significant proinflammatory proteins included in plot:\n")
              print(plot_df$Protein)
              cat("\nTotal count:", nrow(plot_df), "\n")
              
              # ============================================================
              # Bubble plot
              # ============================================================
              bubble_colors <- c(
                "Upregulated (sNPC-MG)"   = "#D73027",
                "Downregulated (sNPC-MG)" = "#4575B4"
              )
              
              n_proteins <- nrow(plot_df)
              plot_height <- max(4, 0.4 * n_proteins + 2)   # auto-scale height by protein count
              
              p_bubble_inflam <- ggplot(plot_df, aes(x = logFC, y = Protein, size = negLogP, color = Category)) +
                geom_point(alpha = 0.85) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
                scale_color_manual(values = bubble_colors) +
                scale_size_continuous(name = expression(bold(paste("-log"[10], "(adj. p-value)"))),
                                      range = c(3, 10)) +
                labs(
                  title = "Significant Proinflammatory Protein Enrichment",
                  subtitle = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = "",
                  color = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(face = "bold", size = 10),
                  legend.position = "right",
                  panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                guides(color = guide_legend(override.aes = list(size = 5)))
              
              print(p_bubble_inflam)
              
              ggsave("ProinflammatoryBubblePlot_AllSignificant_sNPC_vs_sNPCMG.pdf",
                     plot = p_bubble_inflam, width = 8, height = plot_height, dpi = 600)
              ggsave("ProinflammatoryBubblePlot_AllSignificant_sNPC_vs_sNPCMG.png",
                     plot = p_bubble_inflam, width = 8, height = plot_height, dpi = 600)
              
              cat("\nSaved: ProinflammatoryBubblePlot_AllSignificant_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              inflam_df_filtered %>%
                select(Protein, logFC, adj.P.Val, Category) %>%
                arrange(adj.P.Val)
              
              inflam_df_filtered %>%
                dplyr::select(Protein, logFC, adj.P.Val, Category) %>%
                dplyr::arrange(adj.P.Val)
              
              
              
              # ============================================================
              # Bubble Plot: All Matched Proinflammatory Proteins
              # sNPC vs sNPC-MG (significant proteins highlighted)
              # ============================================================
              library(ggplot2)
              library(dplyr)
              
              plot_df <- inflam_df_filtered %>%
                mutate(
                  negLogP = -log10(adj.P.Val),
                  Direction = case_when(
                    logFC > 0 ~ "Up in sNPC-MG",
                    logFC < 0 ~ "Down in sNPC-MG",
                    TRUE ~ "No change"
                  ),
                  Significance = ifelse(adj.P.Val < 0.05, "Significant (adj. p < 0.05)", "Not significant")
                ) %>%
                arrange(logFC)
              
              plot_df$Protein <- factor(plot_df$Protein, levels = plot_df$Protein[order(plot_df$logFC)])
              
              bubble_colors <- c(
                "Up in sNPC-MG"   = "#D73027",
                "Down in sNPC-MG" = "#4575B4",
                "No change"       = "gray60"
              )
              
              p_bubble_inflam <- ggplot(plot_df, aes(x = logFC, y = Protein,
                                                     size = negLogP,
                                                     color = Direction,
                                                     shape = Significance)) +
                geom_point(alpha = 0.85, stroke = 1.1) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
                scale_color_manual(values = bubble_colors) +
                scale_shape_manual(values = c("Significant (adj. p < 0.05)" = 16, "Not significant" = 1)) +
                scale_size_continuous(name = expression(bold(paste("-log"[10], "(adj. p-value)"))),
                                      range = c(3, 12)) +
                labs(
                  title = "Proinflammatory Protein Expression",
                  subtitle = "sNPC vs sNPC-MG",
                  x = expression(bold(paste("log"[2], " Fold Change"))),
                  y = "",
                  color = "",
                  shape = ""
                ) +
                theme_classic(base_size = 13) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "gray30"),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.text.x = element_text(face = "bold", size = 10, color = "black"),
                  axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  axis.ticks = element_line(color = "black", linewidth = 0.6),
                  legend.text = element_text(size = 9),
                  legend.title = element_text(face = "bold", size = 10),
                  legend.position = "right",
                  panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                guides(color = guide_legend(override.aes = list(size = 5)),
                       shape = guide_legend(override.aes = list(size = 5)))
              
              print(p_bubble_inflam)
              
              ggsave("ProinflammatoryBubblePlot_AllMatched_sNPC_vs_sNPCMG.pdf",
                     plot = p_bubble_inflam, width = 8, height = 6, dpi = 600)
              ggsave("ProinflammatoryBubblePlot_AllMatched_sNPC_vs_sNPCMG.png",
                     plot = p_bubble_inflam, width = 8, height = 6, dpi = 600)
              
              cat("Saved: ProinflammatoryBubblePlot_AllMatched_sNPC_vs_sNPCMG.pdf / .png (600 dpi)\n")
              
              
              
              
              ######################## PCA plot ##################
              
              
              library(dplyr)
              library(ggplot2)
              library(ggrepel)
              
              # --- Prepare PCA matrix ---
              pca_matrix <- protein_filtered[complete.cases(protein_filtered), ]
              dim(pca_matrix)
              
              protein_variance <- apply(pca_matrix, 1, var)
              pca_matrix <- pca_matrix[protein_variance > 0, ]
              
              # --- Run PCA (samples as rows) ---
              pca <- prcomp(t(pca_matrix), center = TRUE, scale. = TRUE)
              
              pca_df <- as.data.frame(pca$x)
              pca_df$Sample <- rownames(pca_df)
              
              pca_df <- left_join(pca_df, metadata, by = "Sample")
              
              # --- Variance explained ---
              variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
              variance_explained[1:5]
              
              # --- Plot ---
              p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
                geom_point(size = 5, alpha = 0.85) +
                stat_ellipse(aes(group = Group), type = "norm", linetype = "dashed", linewidth = 0.6) +
                geom_text_repel(size = 4, show.legend = FALSE) +
                labs(
                  title = "PCA of Proteomics Data",
                  x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
                  y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
                  color = "Group"
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  axis.text = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right"
                )
              
              print(p_pca)
              
              ggsave("PCA_Plot_Proteomics.pdf", plot = p_pca, width = 7, height = 6, dpi = 600)
              ggsave("PCA_Plot_Proteomics.png", plot = p_pca, width = 7, height = 6, dpi = 600)
              
              ls()
              
              
              dim(pca_input_matrix)
              head(rownames(pca_input_matrix))
              head(colnames(pca_input_matrix))
              
              str(sample_metadata)   # or meta — check which one has "Sample" and "Group" columns
              
              library(dplyr)
              library(ggplot2)
              library(ggrepel)
              
              pca_matrix <- pca_input_matrix[complete.cases(pca_input_matrix), ]
              dim(pca_matrix)
              
              protein_variance <- apply(pca_matrix, 1, var)
              pca_matrix <- pca_matrix[protein_variance > 0, ]
              
              pca <- prcomp(t(pca_matrix), center = TRUE, scale. = TRUE)
              
              pca_df2 <- as.data.frame(pca$x)
              pca_df2$Sample <- rownames(pca_df2)
              
              pca_df2 <- left_join(pca_df2, sample_metadata, by = "Sample")   # swap sample_metadata -> meta if that's the right one
              
              variance_explained2 <- pca$sdev^2 / sum(pca$sdev^2) * 100
              variance_explained2[1:5]
              
              p_pca <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
                geom_point(size = 5, alpha = 0.85) +
                stat_ellipse(aes(group = Group), type = "norm", linetype = "dashed", linewidth = 0.6) +
                geom_text_repel(size = 4, show.legend = FALSE) +
                labs(
                  title = "PCA of Proteomics Data",
                  x = paste0("PC1 (", round(variance_explained2[1], 1), "%)"),
                  y = paste0("PC2 (", round(variance_explained2[2], 1), "%)"),
                  color = "Group"
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  axis.text = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right"
                )
              
              print(p_pca)
              
              table(pca_df2$Group)
              
              stat_ellipse(aes(group = Group), type = "t", linetype = "dashed", linewidth = 0.6)
              
              
              
              p_pca <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
                geom_point(size = 5, alpha = 0.85) +
                stat_ellipse(aes(group = Group), type = "t", linetype = "dashed", linewidth = 0.6) +
                geom_text_repel(size = 4, show.legend = FALSE) +
                labs(
                  title = "PCA of Proteomics Data",
                  x = paste0("PC1 (", round(variance_explained2[1], 1), "%)"),
                  y = paste0("PC2 (", round(variance_explained2[2], 1), "%)"),
                  color = "Group"
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  axis.text = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right"
                )
              
              print(p_pca)
              
              ggsave("PCA_Plot_Proteomics.pdf", plot = p_pca, width = 7, height = 6, dpi = 600)
              ggsave("PCA_Plot_Proteomics.png", plot = p_pca, width = 7, height = 6, dpi = 600)
              
              
              
              
              ###############elipsed pca
              
              
              # install.packages("ggforce")  # if not already installed
              library(ggforce)
              
              p_pca <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
                geom_mark_ellipse(aes(fill = Group, group = Group), alpha = 0.1, linetype = "dashed",
                                  expand = unit(3, "mm"), show.legend = FALSE) +
                geom_point(size = 5, alpha = 0.85) +
                geom_text_repel(size = 4, show.legend = FALSE) +
                labs(
                  title = "PCA of Proteomics Data",
                  x = paste0("PC1 (", round(variance_explained2[1], 1), "%)"),
                  y = paste0("PC2 (", round(variance_explained2[2], 1), "%)"),
                  color = "Group"
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  axis.text = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right"
                )
              
              print(p_pca)
              
              ggsave("PCA_Plot_Proteomics.pdf", plot = p_pca, width = 7, height = 6, dpi = 600)
              ggsave("PCA_Plot_Proteomics.png", plot = p_pca, width = 7, height = 6, dpi = 600)
              #################
              
              library(ggplot2)
              library(ggrepel)
              library(ggforce)
              
              p_pca <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = Group)) +
                geom_mark_ellipse(aes(fill = Group, group = Group),
                                  alpha = 0.1, linetype = "dashed",
                                  expand = unit(3, "mm"),
                                  show.legend = FALSE) +
                geom_point(size = 5, alpha = 0.85) +
                geom_text_repel(aes(label = Sample), size = 4, show.legend = FALSE) +
                labs(
                  title = "PCA of Proteomics Data",
                  x = paste0("PC1 (", round(variance_explained2[1], 1), "%)"),
                  y = paste0("PC2 (", round(variance_explained2[2], 1), "%)"),
                  color = "Group"
                ) +
                theme_classic(base_size = 14) +
                theme(
                  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  axis.text = element_text(face = "bold", size = 10, color = "black"),
                  axis.line = element_line(color = "black", linewidth = 0.6),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right"
                )
              
              print(p_pca)
              
              ggsave("PCA_Plot_Proteomics.pdf", plot = p_pca, width = 7, height = 6, dpi = 600)
              ggsave("PCA_Plot_Proteomics.png", plot = p_pca, width = 7, height = 6, dpi = 600)
              
              