#Author:Daizee Talukdar

library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(ggpubr)
library(readr)
library (ggtext)
# Load data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Convert location to factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Observed ~ location, data = df)
print(kruskal_res)

# Dunn's test (BH adjusted p-values)
dunn_result <- dunnTest(Observed ~ location, data = df, method = "bh")

# Prepare Dunn results, WITHOUT quotes around asterisks
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# Filter significant comparisons
sig_dunn <- dunn_df %>%
  filter(p_label != "")

# Add y.position manually
max_observed <- max(df$Observed, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_observed + seq(250, by = 70, length.out = n())
  )

# Set factor levels and colors
desired_order <- c("Badkhal", "Ballabhgarh", "Dabua colony", "Parvatiya colony", "Sanjay colony", "BKH", "ESIC")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Badkhal" = "#f25d52",
  "Ballabhgarh" = "#bd900b",
  "Dabua colony" = "#52bf90",
  "Parvatiya colony" = "#cca3ff",
  "Sanjay colony" = "#d962b3",
  "BKH" = "#63a31d",
  "ESIC" = "#11bdd4"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot and save PNG
png(filename = "Observed_fbd_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Observed, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white",
               outlier.size = 1, outlier.fill = "black", outlier.color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Observed", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Enable markdown rendering for p-value labels
    plot.caption = ggtext::element_markdown()
  ) +

  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,   # Must be FALSE for markdown bold
    size = 6         # Adjust size to make bold clearer
  )

dev.off()

##############################################
#Chao1
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(ggtext)  # for markdown element rendering

# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Ensure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Chao1 ~ location, data = df)
print(kruskal_res)

# Dunn's test with BH adjustment
dunn_result <- dunnTest(Chao1 ~ location, data = df, method = "bh")

# Prepare Dunn results without quotes around asterisks
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# Filter significant comparisons only
sig_dunn <- dunn_df %>% filter(p_label != "")

# Add y.position manually for significant comparisons
max_Chao1 <- max(df$Chao1, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_Chao1 + seq(250, by = 70, length.out = n())
  )

# Set factor levels and colors
desired_order <- c("Badkhal", "Ballabhgarh", "Dabua colony", "Parvatiya colony", 
                   "Sanjay colony", "BKH", "ESIC")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Badkhal" = "#f25d52",
  "Ballabhgarh" = "#bd900b",
  "Dabua colony" = "#52bf90",
  "Parvatiya colony" = "#cca3ff",
  "Sanjay colony" = "#d962b3",
  "BKH" = "#63a31d",
  "ESIC" = "#11bdd4"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot and save PNG
png(filename = "Chao1_fbd_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Chao1, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white", 
               outlier.size = 1, outlier.fill = "black", outlier.color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Chao1", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Enable markdown rendering for p-value labels
    plot.caption = ggtext::element_markdown()
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,   # important for markdown bold
    size = 6         # optional size increase for visibility
  )

dev.off()

######################################
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(ggtext)  # for markdown element rendering

# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Make sure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Shannon ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values)
dunn_result <- dunnTest(Shannon ~ location, data = df, method = "bh")

# Prepare Dunn results without quotes around stars
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Filter significant comparisons only
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_Shannon <- max(df$Shannon, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_Shannon + seq(0.7, by = 0.7, length.out = n())
  )

# Set factor levels and colors
desired_order <- c("Badkhal", "Ballabhgarh", "Dabua colony", "Parvatiya colony", "Sanjay colony", "BKH", "ESIC")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Badkhal" = "#f25d52",
  "Ballabhgarh" = "#bd900b",
  "Dabua colony" = "#52bf90",
  "Parvatiya colony" = "#cca3ff",
  "Sanjay colony" = "#d962b3",
  "BKH" = "#63a31d",
  "ESIC" = "#11bdd4"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "Shannon_fbd_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Shannon, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white",
               outlier.size = 1, outlier.fill = "black", outlier.color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Shannon", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Enable markdown rendering for p-value labels
    plot.caption = ggtext::element_markdown()
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,   # <-- important for markdown bold
    size = 6         # adjust label size for visibility
  )

dev.off()

###############################3
#Simpson
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(ggtext)  # for markdown rendering

# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Ensure location is factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Simpson ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values)
dunn_result <- dunnTest(Simpson ~ location, data = df, method = "bh")

# Prepare Dunn results with p_labels without quotes
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Filter significant comparisons only
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_Simpson <- max(df$Simpson, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_Simpson + seq(0.01, by = 0.006, length.out = n())
  )

# Factor levels and colors
desired_order <- c("Badkhal", "Ballabhgarh", "Dabua colony", "Parvatiya colony", "Sanjay colony", "BKH", "ESIC")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Badkhal" = "#f25d52",
  "Ballabhgarh" = "#bd900b",
  "Dabua colony" = "#52bf90",
  "Parvatiya colony" = "#cca3ff",
  "Sanjay colony" = "#d962b3",
  "BKH" = "#63a31d",
  "ESIC" = "#11bdd4"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "Simpson_fbd_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Simpson, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white",
               outlier.size = 1, outlier.fill = "black", outlier.color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Simpson", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = ggtext::element_markdown()
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,   
    size = 6         
  )

dev.off()
