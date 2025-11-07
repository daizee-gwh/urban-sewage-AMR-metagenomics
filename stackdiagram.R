#Author: Daizee Talukdar

library(ggplot2)
library(reshape2)

# convert into relative abundance
library(dplyr)

# Read the file
M <- read.csv("xenobiotics.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Rename first column if needed
names(M)[1] <- "Feature"

# Get only numeric columns (excluding 'Feature')
numeric_cols <- setdiff(names(M), "Feature")

# Convert to numeric (in case some are read as characters)
M[numeric_cols] <- lapply(M[numeric_cols], as.numeric)

# Group by Feature and sum the numeric columns
M_agg <- M %>%
  group_by(Feature) %>%
  summarise_at(vars(numeric_cols), sum, na.rm = TRUE)

# Set Feature as rownames and drop it from data
rownames(M_agg) <- M_agg$Feature
M_agg_1 <- M_agg[, -1]

# Compute relative abundance

# relative abundance as percentages
M_rel <- sweep(M_agg_1, 1, rowSums(M_agg_1), "/") * 100



rel_agg <- data.frame(
  Feature = rownames(M_agg),          # first column = IDs
  M_rel,                              # the numeric abundances
  row.names = NULL                    # drop old row‑name attribute
)


# Optional: round values
M_rel <- round(M_rel, 4)

M_rel <- round(M_rel, 4)

library(tibble)

M_rel_df <- rel_agg %>%
  as.data.frame() %>%              # convert matrix to data frame
  column_to_rownames("Feature")   # add row names as a column called "Feature"

# 1. Load your matrix/data frame
M <- read.csv("xenobiotics.csv", check.names = FALSE)


# Step 2: Keep only numeric columns
M_numeric <- M[, sapply(M, is.numeric)]

# 2. Aggregate rows with the same row names by summing
M_agg <- rowsum(M_numeric, group = rownames(M))

# 3. Compute relative abundance row-wise (each value / row sum)
M1 <- sweep(M_agg, 1, rowSums(M_agg), "/")

dat <- read.table("xenobiotics.csv",
                  header = TRUE,
                  row.names = 1,      # treat first column as row names
                  sep = "\t",         # use "," if CSV
                  check.names = FALSE)

M = read.csv("xenobiotics.csv", header = T)

M1 = sweep(M,1,rowSums(M),"/")
rownames(M_rel) <- M_agg$Feature

plot_df <- rel_agg %>% 
  pivot_longer(
    -Feature,                                # everything except Feature
    names_to  = "Pathway", 
    values_to = "RelAbund"
  )
colnames(rel_agg) <- gsub("\\.", " ", colnames(rel_agg))

colnames(plot_df)


colours <- c("#033CFF",  "#ff0080", "#486935", "#7e508c","#F1A66A","#F26157",
             "#F9ECCC", "#679289", "#33658A","#F6AE2D","#ff5b00","#A54657",
             "#FBB7C2","#ebdb88","#3baa3c","#51C7D9", "#e8f531", "#7519e6",
             "#6e5d41", "#633a59")
plot_df$Feature <- factor(plot_df$Feature,levels=c("Badkhal", "Ballabhgarh", "Dabua colony", "Parvartiya colony", "Sanjay colony", "BKH", "ESIC"))
plot_df$Feature <- factor(plot_df$Feature,levels=c("Haryana", "West Bengal", "Uttar Pradesh", "Uttarakhand", "Assam", "Jharkhand"))
#make the plot!
png(filename = "stack.png",width = 9,height = 8,units = "in",res = 600)
mx = ggplot(plot_df, aes(x = Feature, fill = Pathway, y = RelAbund)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, size = 14, colour = "black", hjust = 1, face= "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "Xenobiotic pathways") + 
  scale_fill_manual(values = colours)

mx
dev.off()
