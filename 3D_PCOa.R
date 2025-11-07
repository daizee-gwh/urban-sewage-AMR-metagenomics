#Author:Daizee Talukdar
# Load necessary libraries
library("phyloseq")
library("vegan")
library("ape")
library("plotly")
library("dplyr")
library("readxl")       
library("dplyr")        
library("tibble") 
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")

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
#normalize
pan_fbd_rel <- transform_sample_counts(pan_fbd, function(x) x / sum(x))
# Step 1: Distance matrix (Bray-Curtis)
dist <- phyloseq::distance(pan_fbd_rel, method = "bray")

# Step 2: PCoA
pcoa_result <- ape::pcoa(dist)

# Calculate percent variance explained by each axis
eig_vals <- pcoa_result$values$Relative_eig
percent_var <- round(100 * eig_vals[1:3], 1)  # First 3 axes

# Extract first 3 coordinates
pcoa_vectors <- as.data.frame(pcoa_result$vectors[, 1:3])
colnames(pcoa_vectors) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_vectors$sample <- rownames(pcoa_vectors)

# Step 3: Convert sample_data to proper data.frame
metadata <- as(sample_data(pan_fbd_rel), "data.frame")
metadata$sample <- rownames(metadata)

# Step 4: Run PERMANOVA
adonis_result <- adonis2(dist ~ Collection_site1, data = metadata, permutations = 999)
print(adonis_result)

# Step 5: Merge ordination with metadata
plot_data <- left_join(pcoa_vectors, metadata, by = "sample")

# Step 6: Define custom hex colors for 6 groups
group_colors <- c(
  "Badkhal" = "firebrick1",
  "Ballabhgarh" = "goldenrod3",
  "Dabua colony" = "#458B74",
  "Parvatiya colony" = "#9A32CD",
  "Sanjay colony" = "#EE1289",
  "BKH" = "#458B00",
  "ESIC" = "deepskyblue2"
)
# Force factor levels to match your colors exactly
plot_data$Collection_site1 <- factor(plot_data$Collection_site1, levels = names(group_colors))

# Now plot with legend title "Collection sites"
p1 <- plot_ly(
  data = plot_data,
  x = ~PCoA1,
  y = ~PCoA2,
  z = ~PCoA3,
  color = ~Collection_site1,
  colors = group_colors,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("<b>PCoA1 (", percent_var[1], "%)</b>")),
      yaxis = list(title = paste0("<b>PCoA2 (", percent_var[2], "%)</b>")),
      zaxis = list(title = paste0("<b>PCoA3 (", percent_var[3], "%)</b>"))
    ),
    title = list(
      text = paste0("<b>PERMANOVA R² = ", round(adonis_result$R2[1], 3),
                    ", p = ", signif(adonis_result$`Pr(>F)`[1], 3), "</b>"),
      xref = "paper",
      x = 0
    ),
    legend = list(
      title = list(text = "<b>Collection sites</b>"),
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )

p1


#############################################3
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
#normalize
pan_wh_rel <- transform_sample_counts(pan_wh, function(x) x / sum(x))
# Step 1: Distance matrix (Bray-Curtis)
dist <- phyloseq::distance(pan_wh_rel, method = "bray")

# Step 2: PCoA
pcoa_result <- ape::pcoa(dist)

# Calculate percent variance explained by each axis
eig_vals <- pcoa_result$values$Relative_eig
percent_var <- round(100 * eig_vals[1:3], 1)  # First 3 axes

# Extract first 3 coordinates
pcoa_vectors <- as.data.frame(pcoa_result$vectors[, 1:3])
colnames(pcoa_vectors) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_vectors$sample <- rownames(pcoa_vectors)

# Step 3: Convert sample_data to proper data.frame
metadata <- as(sample_data(pan_wh_rel), "data.frame")
metadata$sample <- rownames(metadata)

# Step 4: Run PERMANOVA
adonis_result <- adonis2(dist ~ Collection_states, data = metadata, permutations = 999)
print(adonis_result)

# Step 5: Merge ordination with metadata
plot_data <- left_join(pcoa_vectors, metadata, by = "sample")

# Step 6: Define custom hex colors for 6 groups
group_colors <- c(
  "Haryana" = "goldenrod3",
  "Assam" = "#607B8B",
  "Jharkhand" = "#E066FF",
  "West Bengal" = "#EE2C2C",
  "Uttarakhand" = "deepskyblue1",
  "Uttar Pradesh" = "#66CD00"
)

# Step 7: 3D PCoA Plot with percent variance on axes

# Force factor levels to match your colors exactly
plot_data$Collection_site1 <- factor(plot_data$Collection_states, levels = names(group_colors))

# Now plot
p2 <- plot_ly(
  data = plot_data,
  x = ~PCoA1,
  y = ~PCoA2,
  z = ~PCoA3,
  color = ~Collection_states,
  colors = group_colors,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("<b>PCoA1 (", percent_var[1], "%)</b>")),
      yaxis = list(title = paste0("<b>PCoA2 (", percent_var[2], "%)</b>")),
      zaxis = list(title = paste0("<b>PCoA3 (", percent_var[3], "%)</b>"))
    ),
    title = list(
      text = paste0("<b>PERMANOVA R² = ", round(adonis_result$R2[1], 3),
                    ", p = ", signif(adonis_result$`Pr(>F)`[1], 3), "</b>"),
      xref = "paper",
      x = 0
    ),
    legend = list(
      title = list(text = "<b>Collection states</b>"),  # <-- Legend title added here
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )

p2


############################################################
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

#normalize
pan_m_rel <- transform_sample_counts(pan_m, function(x) x / sum(x))

# Step 1: Distance matrix (Bray-Curtis)
dist <- phyloseq::distance(pan_m_rel, method = "bray")

# Step 2: PCoA
pcoa_result <- ape::pcoa(dist)

# Calculate percent variance explained by each axis
eig_vals <- pcoa_result$values$Relative_eig
percent_var <- round(100 * eig_vals[1:3], 1)  # First 3 axes

# Extract first 3 coordinates
pcoa_vectors <- as.data.frame(pcoa_result$vectors[, 1:3])
colnames(pcoa_vectors) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_vectors$sample <- rownames(pcoa_vectors)

# Step 3: Convert sample_data to proper data.frame
metadata <- as(sample_data(pan_m_rel), "data.frame")
metadata$sample <- rownames(metadata)

# Step 4: Run PERMANOVA
adonis_result <- adonis2(dist ~ Months, data = metadata, permutations = 999)
print(adonis_result)

# Step 5: Merge ordination with metadata
plot_data <- left_join(pcoa_vectors, metadata, by = "sample")

# Step 6: Define custom hex colors for 6 groups
group_colors <- c(
  "June" = "firebrick1",
  "July" = "goldenrod3",
  "August" = "olivedrab3",
  "September" = "#43CD80",
  "October" = "deepskyblue2",
  "November" = "#9A32CD",
  "December" = "#FF34B3"
)

# Force factor levels to match your colors exactly
plot_data$Months <- factor(plot_data$Months, levels = names(group_colors))

# Now plot
p3 <- plot_ly(
  data = plot_data,
  x = ~PCoA1,
  y = ~PCoA2,
  z = ~PCoA3,
  color = ~Months,
  colors = group_colors,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("<b>PCoA1 (", percent_var[1], "%)</b>")),
      yaxis = list(title = paste0("<b>PCoA2 (", percent_var[2], "%)</b>")),
      zaxis = list(title = paste0("<b>PCoA3 (", percent_var[3], "%)</b>"))
    ),
    title = list(
      text = paste0("<b>PERMANOVA R² = ", round(adonis_result$R2[1], 3),
                    ", p = ", signif(adonis_result$`Pr(>F)`[1], 3), "</b>"),
      xref = "paper",
      x = 0
    ),
    legend = list(
      title = list(text = "<b>Months</b>"),  # Legend title added here
      font = list(
        family = "Arial Black, Arial, sans-serif",
        size = 12,
        color = "black"
      )
    )
  )

p3
