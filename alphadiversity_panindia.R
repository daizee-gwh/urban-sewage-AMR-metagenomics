#Author:Daizee Talukdar

library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(ggpubr)
library(readr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(ggtext)  # for markdown rendering

# Load your data
df <- read_csv("alphadiversity_wh.csv")

# Make sure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Observed ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values)
dunn_result <- dunnTest(Observed ~ location, data = df, method = "bh")

# Prepare Dunn results
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

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_observed <- max(df$Observed, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_observed + seq(170, by = 76, length.out = n())
  )

# Factor levels and colors
desired_order <- c("Haryana", "West Bengal", "Uttar Pradesh", "Uttarakhand", "Assam", "Jharkhand")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Haryana" = "#bd900b",
  "West Bengal" = "#f25d52",
  "Uttar Pradesh" = "#63a31d",
  "Uttarakhand" = "#11bdd4",
  "Assam" = "#7e9cab",
  "Jharkhand" = "#d962b3"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "Observed_wh_wg5.png", width = 7, height = 6, units = "in", res = 600)

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
    parse = FALSE,   # Important for showing stars correctly
    size = 6         # Adjust size as needed
  )

dev.off()

############################ shannon
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(ggpubr)
library(readr)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(ggtext)  # for markdown rendering if needed

# Load your data
df <- read_csv("alphadiversity_wh.csv")

# Make sure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Shannon ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values)
dunn_result <- dunnTest(Shannon ~ location, data = df, method = "bh")

# Prepare Dunn results
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

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_shannon <- max(df$Shannon, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_shannon + seq(0.5, by = 0.5, length.out = n())  # Adjust vertical spacing here as needed
  )

# Factor levels and colors
desired_order <- c("Haryana", "West Bengal", "Uttar Pradesh", "Uttarakhand", "Assam", "Jharkhand")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Haryana" = "#bd900b", "West Bengal" = "#f25d52", "Uttar Pradesh" = "#63a31d",
  "Uttarakhand" = "#11bdd4", "Assam" = "#7e9cab", "Jharkhand" = "#d962b3"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "shannon_wh_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Shannon, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white", outlier.size = 1,
               outlier.fill = "black", outlier.color = "black") +
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
    plot.caption = ggtext::element_markdown()  # if you want markdown rendering
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,  # Show stars as is, no parsing needed
    size = 6        # Adjust size if you want bigger stars
  )

dev.off()

#################################Chao1
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)

# Load your data
df <- read_csv("alphadiversity_wh.csv")

# Make sure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test for Chao1
kruskal_res <- kruskal.test(Chao1 ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values) for Chao1
dunn_result <- dunnTest(Chao1 ~ location, data = df, method = "bh")

# Prepare Dunn results
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

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>% filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_chao1 <- max(df$Chao1, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_chao1 + seq(170, by = 76, length.out = n())  # Keep your spacing
  )

# Factor levels and colors
desired_order <- c("Haryana", "West Bengal", "Uttar Pradesh", "Uttarakhand", "Assam", "Jharkhand")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Haryana" = "#bd900b", "West Bengal" = "#f25d52", "Uttar Pradesh" = "#63a31d",
  "Uttarakhand" = "#11bdd4", "Assam" = "#7e9cab", "Jharkhand" = "#d962b3"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "Chao1_wh_dunn_wg5.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Chao1, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white", outlier.size = 1,
               outlier.fill = "black", outlier.color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Chao1", x = "") +
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
    panel.grid.minor = element_blank()
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    parse = FALSE,   # Ensure stars display correctly
    size = 6         # Optional: bigger labels
  )

dev.off()

################################## Simpson
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FSA)

# Load your data
df <- read_csv("alphadiversity_wh.csv")

# Make sure location is a factor
df$location <- as.factor(df$location)

# Kruskal-Wallis test for Simpson
kruskal_res <- kruskal.test(Simpson ~ location, data = df)
print(kruskal_res)

# Dunn's test (adjusted p-values) for Simpson
dunn_result <- dunnTest(Simpson ~ location, data = df, method = "bh")

# Prepare Dunn results
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

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>% filter(p_label != "ns")

# Add y.position manually for significant comparisons
max_simpson <- max(df$Simpson, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = max_simpson + seq(0.02, by = 0.006, length.out = n())  # Keep your spacing
  )

# Factor levels and colors
desired_order <- c("Haryana", "West Bengal", "Uttar Pradesh", "Uttarakhand", "Assam", "Jharkhand")
df$location <- factor(df$location, levels = desired_order)

custom_colors <- c(
  "Haryana" = "#bd900b", "West Bengal" = "#f25d52", "Uttar Pradesh" = "#63a31d",
  "Uttarakhand" = "#11bdd4", "Assam" = "#7e9cab", "Jharkhand" = "#d962b3"
)
custom_colors <- custom_colors[levels(df$location)]

# Plot
png(filename = "Simpson_wh_dunn_wg5.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = location, y = Simpson, fill = location)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = 21, fill = "white", outlier.size = 1,
               outlier.fill = "black", outlier.color = "black") +
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
    panel.grid.minor = element_blank()
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
