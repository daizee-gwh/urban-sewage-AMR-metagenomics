#Author:Daizee Talukdar
#beta diversity calculation for community and hospital in accordance with month 

#metadata, dada2 output sheet, taxonomy
# create a phyloseq object
library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(tibble) 
library(ggpubr)
library(Matrix)
library(reshape2)
library(vegan)
library(plotly)
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
#normalize
site_fbd1 <- transform_sample_counts(site_fbd, function(x) x / sum(x))
# Bray-Curtis distance on raw counts (no transformation)
bray_dist <- distance(site_fbd1, method = "bray")

# Run 3D NMDS
set.seed(123)
nmds_3d <- metaMDS(bray_dist, k = 3, trymax = 100)

# Extract stress value
stress_value <- nmds_3d$stress

# Extract NMDS scores (3D)
nmds_scores <- as.data.frame(scores(nmds_3d))
nmds_scores$SampleID <- rownames(nmds_scores)

# Prepare metadata dataframe from phyloseq
metadata_df <- as(sample_data(site_fbd1), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# Merge NMDS scores with metadata
nmds_plot_data <- merge(nmds_scores, metadata_df, by = "SampleID")

# PERMANOVA test using adonis2 (vegan)
adonis_results <- adonis2(bray_dist ~ Settings * Months, data = metadata_df)
# Extract values for annotation
R2 <- round(adonis_results$R2[1], 4)
Fvalue <- round(adonis_results$F[1], 4)
pvalue <- adonis_results$`Pr(>F)`[1]

# Custom colors for Settings
settings_colors <- c("Community" = "#d24b4b", "Hospital" = "#74DDD9")

# Your unique months in the dataset
unique_months <- unique(nmds_plot_data$Months)

# Define month shapes only for months present in your data
month_shapes <- c(
  "M2_July" = "circle", 
  "M3_August" = "square", 
  "M4_September" = "diamond",
  "M5_October" = "diamond-open", 
  "M6_November" ="circle-open",
  "M7_December" = "square-open"
)

# Subset month_shapes to only those present in your data (just in case)
month_shapes <- month_shapes[names(month_shapes) %in% unique_months]

# Make sure Months is a factor with levels exactly matching the shape keys
nmds_plot_data$Months <- factor(nmds_plot_data$Months, levels = names(month_shapes))

# ----------------------------
# Font Settings
# ----------------------------

# Clean font for title and axis titles
title_font <- list(family = "Arial", size = 14, color = "black")
axis_title_font <- list(family = "Arial", size = 14, color = "black")

# Bolder font for legend title and text
legend_font <- list(
  family = "Arial Black, Arial, sans-serif",
  size = 14,
  color = "black"
)

# ----------------------------
# Compose NMDS plot title
# ----------------------------

full_title <- paste0(
  "<b>3D NMDS Plot (Bray-Curtis)<br>",
  "Stress: ", round(stress_value, 4),
  " | PERMANOVA R²: ", R2,
  " F: ", Fvalue,
  " p: ", signif(pvalue, 3), "</b>"
)

# ----------------------------
# Plot
# ----------------------------

p <- plot_ly(
  data = nmds_plot_data,
  x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
  type = 'scatter3d',
  mode = 'markers',
  color = ~Settings,
  colors = settings_colors,
  symbol = ~Months,
  symbols = month_shapes,
  marker = list(size = 5)
) %>%
  layout(
    title = list(
      text = full_title,
      font = title_font
    ),
    scene = list(
      xaxis = list(title = list(text = "<b>NMDS1</b>", font = axis_title_font)),
      yaxis = list(title = list(text = "<b>NMDS2</b>", font = axis_title_font)),
      zaxis = list(title = list(text = "<b>NMDS3</b>", font = axis_title_font))
    ),
    legend = list(
      title = list(
        text = "Settings & Months",
        font = legend_font  
      ),
      font = legend_font     
    )
  )

# Show the plot
p
