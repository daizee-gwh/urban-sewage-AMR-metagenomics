#Author: Daizee Talukdar
# Load Required Libraries 
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(igraph)
library(ggraph)


setwd("~july")  

#  Load and Normalize Data 
# Bacterial genes
bac_data <- read.csv("final.csv", row.names = 1, stringsAsFactors = FALSE, check.names = F)
bac_rel <- sweep(as.matrix(bac_data), 1, rowSums(bac_data), FUN = "/")
bac_rel <- t(bac_rel)  # genes x samples

# LCMS (antibiotic/metabolite data)
lcms_data <- read.csv("lcms.csv", row.names = 1,check.names = F)
lcms_rel <- sweep(as.matrix(t(lcms_data)), 1, rowSums(t(lcms_data)), FUN = "/")

# Match Sample Columns 
common_samples <- intersect(colnames(bac_rel), colnames(lcms_rel))
bac_rel <- bac_rel[, common_samples]
lcms_rel <- lcms_rel[, common_samples]

# Clean column names (remove trailing dots) 
colnames(bac_rel) <- gsub("\\.$", "", colnames(bac_rel))
colnames(lcms_rel) <- gsub("\\.$", "", colnames(lcms_rel))
rownames(bac_rel) <- gsub("\\.$", "", rownames(bac_rel))
rownames(lcms_rel) <- gsub("\\.$", "", rownames(lcms_rel))

# Spearman Correlation 
cor_matrix <- cor(t(bac_rel), t(lcms_rel), method = "spearman", use = "pairwise.complete.obs")

# Compute p-values 
get_pval <- function(x, y) {
  tryCatch(cor.test(x, y, method = "spearman")$p.value, error = function(e) NA)
}
pval_mat <- outer(rownames(cor_matrix), colnames(cor_matrix), Vectorize(function(i, j) {
  get_pval(bac_rel[i, ], lcms_rel[j, ])
}))
rownames(pval_mat) <- rownames(cor_matrix)
colnames(pval_mat) <- colnames(cor_matrix)

# Save Outputs 
write.csv(cor_matrix, "gene_antibiotic_spearman_correlation_matrix.csv", quote = FALSE)
write.csv(pval_mat, "gene_antibiotic_spearman_pvalues.csv", quote = FALSE)

# Flatten and Save Table 
cor_table <- as.data.frame(as.table(cor_matrix))
colnames(cor_table) <- c("Gene", "Antibiotic", "Spearman_r")

pval_table <- as.data.frame(as.table(pval_mat))
colnames(pval_table) <- c("Gene", "Antibiotic", "p_value")

full_table <- left_join(cor_table, pval_table, by = c("Gene", "Antibiotic"))
write.csv(full_table, "gene_antibiotic_spearman_correlation_table.csv", row.names = FALSE)

sig_table <- full_table %>% filter(abs(Spearman_r) > 0.2 & p_value < 0.05)
write.csv(sig_table, "gene_antibiotic_spearman_table_sig.csv", row.names = FALSE)

# Prepare Heatmap 
cor_clean <- cor_matrix
cor_clean[is.na(cor_clean) | is.nan(cor_clean) | is.infinite(cor_clean)] <- 0

# Custom diverging color palette (red to white to blue)
custom_col <- colorRamp2(
  breaks = seq(-1, 1, length.out = 11),
  colors = colorRampPalette(c("red3", "white", "navy"))(11)
)

# Significant correlation mask
sig_mask <- abs(cor_matrix) > 0.2 & pval_mat < 0.05
asterisk_mat <- ifelse(sig_mask, "*", "")


library(RColorBrewer)
library(circlize)

# Use 11 evenly spaced values from −1 to +1
breaks <- seq(-1, 1, length.out = 11)
colors <- rev(brewer.pal(11, "RdBu"))  # Reversed to make blue = positive, red = negative
cor_col <- colorRamp2(breaks, colors)


#  ComplexHeatmap Plot 
png("gene_antibiotic_spearman_heatmap2.png", width = 2400, height = 2400, res = 300)

Heatmap(cor_clean,
        name = "Spearman r",
        col = cor_col,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontface = "bold.italic"),
        column_names_gp = gpar(fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (!is.na(asterisk_mat[i, j]) && asterisk_mat[i, j] == "*") {
            grid.text("*", x = x, y = y, gp = gpar(fontface = "bold", fontsize = 10))
          }
        })

dev.off()

# Optional: Correlation Network 
graph_edges <- sig_table %>% select(Gene, Antibiotic, Spearman_r) %>%
  rename(from = Gene, to = Antibiotic, weight = Spearman_r)

graph <- graph_from_data_frame(graph_edges, directed = FALSE)

png("gene_antibiotic_spearman_network2.png", width = 2400, height = 2400, res = 300)
ggraph(graph) +
  geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight), color = weight)) +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("brown2", "blue3")) +
  geom_node_point(color = "gold", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold") +
  theme_graph()
dev.off()


