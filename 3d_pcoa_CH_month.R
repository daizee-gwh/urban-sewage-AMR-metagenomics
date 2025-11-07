#Author:Daizee Talukdar
#beta diversity calculation for community and hospital in accordance with month 

#metadata, dada2 output sheet, taxonomy
setwd("/monthwise")
# create a phyloseq object
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble") 
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")
library("ape")
library("plotly")
# Needed for converting column to row names
otu_mat<- read_excel("otu.xlsx", sheet = "otu")
tax_mat<- read_excel("taxonomy.xlsx", sheet = "taxonomy")
samples_df <- read_excel("metadata.xlsx", sheet = "metadata")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
site_fbd <- phyloseq(OTU, TAX, samples)
site_fbd

# # CONVERT TO RELATIVE ABUNDANCE
site_fbd_rel <- transform_sample_counts(site_fbd, function(x) x / sum(x))
# 1. Calculate Bray-Curtis distance matrix
bray_dist <- distance(site_fbd_rel , method = "bray")

# 2. Run PCoA
pcoa_res <- ape::pcoa(bray_dist)

# Calculate % variance explained by first 3 axes
var_exp <- round(100 * pcoa_res$values$Relative_eig[1:3], 1)

# Extract first 3 PCoA axes (coordinates)
pcoa_coords <- as.data.frame(pcoa_res$vectors[, 1:3])
colnames(pcoa_coords) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_coords$SampleID <- rownames(pcoa_coords)

# 3. Prepare metadata dataframe
metadata_df <- as(sample_data(site_fbd_rel), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# 4. Merge PCoA coords and metadata
pcoa_plot_data <- merge(pcoa_coords, metadata_df, by = "SampleID")

# 5. PERMANOVA test (adonis2)
adonis_results <- adonis2(bray_dist ~ Settings * Months, data = metadata_df)
R2 <- round(adonis_results$R2[1], 4)
Fvalue <- round(adonis_results$F[1], 4)
pvalue <- adonis_results$`Pr(>F)`[1]

# 6. Define colors for Settings
settings_colors <- c("Community" = "#d24b4b", "Hospital" = "#74DDD9")

# 7. Define shapes for Months (adjust to your months)
month_shapes <- c("M2_July" = "circle", 
                  "M3_August" = "square", 
                  "M4_September" = "diamond",
                  "M5_October" = "diamond-open", 
                  "M6_November" ="circle-open",
  "M7_December" = "square-open")

# Filter month_shapes to only months present in your data
month_shapes <- month_shapes[names(month_shapes) %in% unique(pcoa_plot_data$Months)]

# 8. Create 3D PCoA plot with bold text and PERMANOVA info
p4 <- plot_ly(data = pcoa_plot_data,
             x = ~PCoA1, y = ~PCoA2, z = ~PCoA3,
             type = 'scatter3d',
             mode = 'markers',
             color = ~Settings,
             colors = settings_colors,
             symbol = ~Months,
             symbols = month_shapes,
             marker = list(size = 5)) %>%
  layout(
    title = list(
      text = paste0("3D PCoA Plot (Bray-Curtis)<br>",
                    "PERMANOVA R²: ", R2, " | F: ", Fvalue, " | p = ", signif(pvalue, 3)),
      font = list(family = "Arial Black", size = 16)
    ),
    scene = list(
      xaxis = list(title = list(text = paste0("PCoA1 (", var_exp[1], "%)"),
                                font = list(family = "Arial Black", size = 14))),
      yaxis = list(title = list(text = paste0("PCoA2 (", var_exp[2], "%)"),
                                font = list(family = "Arial Black", size = 14))),
      zaxis = list(title = list(text = paste0("PCoA3 (", var_exp[3], "%)"),
                                font = list(family = "Arial Black", size = 14)))
    ),
    legend = list(
      title = list(text = "<b>Settings & Months</b>"),
      font = list(size = 12)
    )
  )

# 9. Show the plot
p4
