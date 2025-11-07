#Author:Daizee Talukdar
# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(plotly)
library(dplyr)
library(readxl)      
library(dplyr)        
library(tibble) 
library(ggpubr)
library(Matrix)
library(reshape2)
library(vegan)
setwd("")
#make the phyloseq object for only faridabad 
# Needed for converting column to row names
otu_mat<- read_excel("otu.xlsx", sheet = "otu")
tax_mat<- read_excel("taxonomy.xlsx", sheet = "taxonomy")
samples_df <- read_excel("metadata.xlsx", sheet = "faridabad")
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
pan_fbd <- phyloseq(OTU, TAX, samples)
pan_fbd
# CONVERT TO RELATIVE ABUNDANCE
pan_fbd_rel <- transform_sample_counts(pan_fbd, function(x) x / sum(x))
# Step 1: Calculate Bray-Curtis distance
# -------------------------------
bray_dist <- phyloseq::distance(pan_fbd_rel, method = "bray")

# -------------------------------
# Step 2: Perform NMDS (3D)
# -------------------------------
set.seed(123)  # for reproducibility
nmds <- metaMDS(bray_dist, k = 3, trymax = 100)

# -------------------------------
# Step 3: Extract NMDS coordinates + merge with sample data
# -------------------------------
nmds_coords <- as.data.frame(nmds$points)
nmds_coords$sample <- rownames(nmds_coords)

sample_df <- data.frame(sample_data(pan_fbd_rel))
sample_df$sample <- rownames(sample_df)

plot_data <- merge(nmds_coords, sample_df, by = "sample")

# -------------------------------
# Step 4: Define hex colors for groups
# -------------------------------
grouping_var <- "Collection_site1"

group_levels <- unique(plot_data[[grouping_var]])
group_colors <- c("Badkhal" = "firebrick1",
                  "Ballabhgarh" = "goldenrod3",
                  "Dabua colony" = "#458B74",
                  "Parvatiya colony" = "#9A32CD",
                  "Sanajy colony" = "#EE1289",
                  "BKH" = "#458B00",
                  "ESIC" = "deepskyblue2"
)[seq_along(group_levels)]
names(group_colors) <- group_levels

plot_data[[grouping_var]] <- factor(plot_data[[grouping_var]], levels = group_levels)

# -------------------------------
# Step 5: Plot 3D NMDS with Plotly (with stress and R² in title)
# -------------------------------

# Extract stress value
stress_val <- round(nmds$stress, 4)

# Run PERMANOVA and extract R²
adonis_result <- adonis2(bray_dist ~ get(grouping_var), data = sample_df, permutations = 999)
r2 <- round(adonis_result$R2[1], 3)
#plot 
fig <- plot_ly(
  data = plot_data,
  x = ~MDS1, y = ~MDS2, z = ~MDS3,
  type = 'scatter3d',
  mode = 'markers',
  color = ~get(grouping_var),
  colors = group_colors,
  marker = list(size = 5),
  text = ~sample
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(text = "<b>NMDS1</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      yaxis = list(
        title = list(text = "<b>NMDS2</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      zaxis = list(
        title = list(text = "<b>NMDS3</b>", font = list(family = "Arial", size = 14, color = "black"))
      )
    ),
    title = list(
      text = paste0("<b>3D NMDS | Stress = ", round(stress_val, 3),
                    " | PERMANOVA R² = ", round(r2, 3), "</b>"),
      font = list(family = "Arial", size = 16)
    ),
    legend = list(
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )
# Display plot
fig
# -------------------------------
# Step 6: Print PERMANOVA Results
# -------------------------------
print(adonis_result)

#############################
#without hospital
# removing hospitals
otu_mat<- read_excel("otu.xlsx", sheet = "otu")
tax_mat<- read_excel("taxonomy.xlsx", sheet = "taxonomy")
samples_df <- read_excel("Metadata.xlsx", sheet = "Metadata")
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
pan_wh <- phyloseq(OTU, TAX, samples)
pan_wh
# CONVERT TO RELATIVE ABUNDANCE
pan_wh_rel <- transform_sample_counts(pan_wh, function(x) x / sum(x))
# -------------------------------
# Step 1: Calculate Bray-Curtis distance
# -------------------------------
bray_dist <- phyloseq::distance(pan_wh_rel, method = "bray")

# -------------------------------
# Step 2: Perform NMDS (3D)
# -------------------------------
set.seed(123)  # for reproducibility
nmds <- metaMDS(bray_dist, k = 3, trymax = 100)

# -------------------------------
# Step 3: Extract NMDS coordinates + merge with sample data
# -------------------------------
nmds_coords <- as.data.frame(nmds$points)
nmds_coords$sample <- rownames(nmds_coords)

sample_df <- data.frame(sample_data(pan_wh_rel))
sample_df$sample <- rownames(sample_df)

plot_data <- merge(nmds_coords, sample_df, by = "sample")

# -------------------------------
# Step 4: Define hex colors for groups
# -------------------------------
grouping_var <- "Collection_states"

# Sort states alphabetically
group_levels <- sort(unique(plot_data[[grouping_var]]))  # Alphabetical order

# Define matching colors
group_colors <- c("Assam" = "#607B8B",
                  "Haryana" = "goldenrod3",
                  "Jharkhand" = "#E066FF",
                  "Uttarakhand" = "deepskyblue1",
                  "Uttar Pradesh" = "#66CD00",
                  "West Bengal" = "#EE2C2C"
)[group_levels]  # This ensures the order matches sorted levels

# Set factor levels in alphabetical order
plot_data[[grouping_var]] <- factor(plot_data[[grouping_var]], levels = group_levels)


# -------------------------------
# Step 5: Plot 3D NMDS with Plotly (with stress and R² in title)
# -------------------------------

# Extract stress value
stress_val <- round(nmds$stress, 4)

# Run PERMANOVA and extract R²
adonis_result <- adonis2(bray_dist ~ get(grouping_var), data = sample_df, permutations = 999)
r2 <- round(adonis_result$R2[1], 3)

#plot 
fig2 <- plot_ly(
  data = plot_data,
  x = ~MDS1, y = ~MDS2, z = ~MDS3,
  type = 'scatter3d',
  mode = 'markers',
  color = ~get(grouping_var),
  colors = group_colors,
  marker = list(size = 5),
  text = ~sample
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(text = "<b>NMDS1</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      yaxis = list(
        title = list(text = "<b>NMDS2</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      zaxis = list(
        title = list(text = "<b>NMDS3</b>", font = list(family = "Arial", size = 14, color = "black"))
      )
    ),
    title = list(
      text = paste0("<b>3D NMDS | Stress = ", round(stress_val, 3),
                    " | PERMANOVA R² = ", round(r2, 3), "</b>"),
      font = list(family = "Arial", size = 16)
    ),
    legend = list(
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )
# Display plot
fig2
# -------------------------------
# Step 6: Print PERMANOVA Results
# -------------------------------
print(adonis_result)
##################################
# for months june to december
#make the phyloseq object for only faridabad 
# Needed for converting column to row names
otu_mat<- read_excel("otu.xlsx", sheet = "otu")
tax_mat<- read_excel("taxonomy.xlsx", sheet = "taxonomy")
samples_df <- read_excel("metadata.xlsx", sheet = "faridabad")
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
pan_m <- phyloseq(OTU, TAX, samples)
pan_m
# # CONVERT TO RELATIVE ABUNDANCE
pan_m_rel <- transform_sample_counts(pan_m, function(x) x / sum(x))
# Load Required Libraries
# -------------------------------
library(phyloseq)
library(vegan)
library(plotly)
library(dplyr)
# ------------------------------------------
# Step 1: Calculate Bray-Curtis Distance Matrix
# ------------------------------------------
bray_dist <- phyloseq::distance(pan_m_rel, method = "bray")

# ------------------------------------------
# Step 2: Perform NMDS (3D)
# ------------------------------------------
set.seed(123)
nmds <- metaMDS(bray_dist, k = 3, trymax = 100)

# ------------------------------------------
# Step 3: Extract NMDS Coordinates & Merge with Metadata
# ------------------------------------------
nmds_coords <- as.data.frame(nmds$points)
nmds_coords$sample <- rownames(nmds_coords)

sample_df <- data.frame(sample_data(pan_m_rel))
sample_df$sample <- rownames(sample_df)

plot_data <- merge(nmds_coords, sample_df, by = "sample")

# ------------------------------------------
# Step 4: Define Grouping Variable and Custom Colors
# ------------------------------------------
grouping_var <- "Months"

group_colors <- c(
  "June" = "firebrick1",
  "July" = "goldenrod3",
  "August" = "olivedrab3",
  "September" = "#43CD80",
  "October" = "deepskyblue2",
  "November" = "#9A32CD",
  "December" = "#FF34B3"
)

# Ensure factor levels match color order
plot_data[[grouping_var]] <- factor(plot_data[[grouping_var]], levels = names(group_colors))

# ------------------------------------------
# Step 5: Run PERMANOVA
# ------------------------------------------
adonis_result <- adonis2(bray_dist ~ Months, data = sample_df, permutations = 999)
r2 <- round(adonis_result$R2[1], 3)
pval <- signif(adonis_result$`Pr(>F)`[1], 3)
stress_val <- round(nmds$stress, 4)

# ------------------------------------------
# Step 6: Plot 3D NMDS with Plotly
# ------------------------------------------
fig3 <- plot_ly(
  data = plot_data,
  x = ~MDS1, y = ~MDS2, z = ~MDS3,
  type = 'scatter3d',
  mode = 'markers',
  color = ~get(grouping_var),
  colors = group_colors,
  marker = list(size = 5),
  text = ~sample
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(text = "<b>NMDS1</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      yaxis = list(
        title = list(text = "<b>NMDS2</b>", font = list(family = "Arial", size = 14, color = "black"))
      ),
      zaxis = list(
        title = list(text = "<b>NMDS3</b>", font = list(family = "Arial", size = 14, color = "black"))
      )
    ),
    title = list(
      text = paste0("<b>3D NMDS | Stress = ", stress_val,
                    " | PERMANOVA R² = ", r2, ", p = ", pval, "</b>"),
      font = list(family = "Arial", size = 16)
    ),
    legend = list(
      orientation = "v",
      x = 1,
      xanchor = "left",
      y = 1,
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )

# ------------------------------------------
# Step 7: Display Plot
# ------------------------------------------
fig3

# ------------------------------------------
# Step 8: Print PERMANOVA Table
# ------------------------------------------
print(adonis_result)
