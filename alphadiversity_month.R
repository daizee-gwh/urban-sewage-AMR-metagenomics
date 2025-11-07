#Author:Daizee Talukdar
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(ggpubr)
library(readr)

# Load your data
setwd("")

# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Make sure Months is a factor with desired order
desired_order <- c("June", "July", "August", "September", "October", "November", "December")
df$Months <- factor(df$Months, levels = desired_order)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Observed ~ Months, data = df)
print(kruskal_res)

# Dunn's test with BH correction
dunn_result <- dunnTest(Observed ~ Months, data = df, method = "bh")

# Prepare Dunn results
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = ifelse(p_value < 0.001, "***",
                     ifelse(p_value < 0.01, "**",
                            ifelse(p_value < 0.05, "*", "ns")))
  )

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position for p-value lines
max_observed <- max(df$Observed, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = seq(from = max_observed * 1.05, by = max_observed * 0.05, length.out = n())
  )

# Custom colors for Months
custom_colors <- c(
  "June" = "firebrick1",
  "July" = "goldenrod3",
  "August" = "olivedrab3",
  "September" = "#43CD80",
  "October" = "deepskyblue2",
  "November" = "#9A32CD",
  "December" = "#FF34B3"
)
custom_colors <- custom_colors[levels(df$Months)]

# Plot only filled boxplots with custom colors

png(filename = "Observed_months_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = Months, y = Observed, fill = Months)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 21, outlier.fill = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Observed", x = "") +
  theme_minimal() + 
  stat_summary(fun.data = mean_sd, geom = "errorbar", width= 0.2, colour = "black") +

  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"), 
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )

dev.off()

###############################Shannon
# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Make sure Months is a factor with desired order
desired_order <- c("June", "July", "August", "September", "October", "November", "December")
df$Months <- factor(df$Months, levels = desired_order)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Shannon ~ Months, data = df)
print(kruskal_res)

# Dunn's test with BH correction
dunn_result <- dunnTest(Shannon ~ Months, data = df, method = "bh")

# Prepare Dunn results
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = ifelse(p_value < 0.001, "***",
                     ifelse(p_value < 0.01, "**",
                            ifelse(p_value < 0.05, "*", "ns")))
  )

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position for p-value lines
max_shannon <- max(df$Shannon, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = seq(from = max_shannon * 1.05, by = max_shannon * 0.05, length.out = n())
  )

# Custom colors for Months
custom_colors <- c(
  "June" = "firebrick1",
  "July" = "goldenrod3",
  "August" = "olivedrab3",
  "September" = "#43CD80",
  "October" = "deepskyblue2",
  "November" = "#9A32CD",
  "December" = "#FF34B3"
)
custom_colors <- custom_colors[levels(df$Months)]

# Plot only filled boxplots with custom colors

png(filename = "Shannon_months_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = Months, y = Shannon, fill = Months)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 21, outlier.fill = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = "Shannon", x = "") +
  theme_minimal() + 
  stat_summary(fun.data = mean_sd, geom = "errorbar", width= 0.2, colour = "black") +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"), 
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )

dev.off()

########################################sIMPSON

# Load your data
df <- read_csv("Diversity_metadatasamples_updated.csv")

# Make sure Months is a factor with desired order
desired_order <- c("June", "July", "August", "September", "October", "November", "December")
df$Months <- factor(df$Months, levels = desired_order)

# Kruskal-Wallis test
kruskal_res <- kruskal.test(Simpson ~ Months, data = df)
print(kruskal_res)

# Dunn's test with BH correction
dunn_result <- dunnTest(Simpson ~ Months, data = df, method = "bh")

# Prepare Dunn results
dunn_df <- dunn_result$res %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  rename(p_value = P.adj) %>%
  mutate(
    p_label = ifelse(p_value < 0.001, "***",
                     ifelse(p_value < 0.01, "**",
                            ifelse(p_value < 0.05, "*", "ns")))
  )

# Filter only significant comparisons (p < 0.05)
sig_dunn <- dunn_df %>%
  filter(p_label != "ns")

# Add y.position for p-value lines
max_simpson <- max(df$Simpson, na.rm = TRUE)
sig_dunn <- sig_dunn %>%
  mutate(
    y.position = seq(from = max_simpson * 1.01, by = 0.002, length.out = n())  # space between p-value lines
  )

# Custom colors for Months
custom_colors <- c(
  "June" = "firebrick1",
  "July" = "goldenrod3",
  "August" = "olivedrab3",
  "September" = "#43CD80",
  "October" = "deepskyblue2",
  "November" = "#9A32CD",
  "December" = "#FF34B3"
)
custom_colors <- custom_colors[levels(df$Months)]

# Plot only filled boxplots with custom colors
png(filename = "Simpson_month_dunn.png", width = 7, height = 6, units = "in", res = 600)

ggplot(df, aes(x = Months, y = Simpson, fill = Months)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 21, outlier.color = "black", outlier.fill = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(
    limits = c(0.98, max(sig_dunn$y.position) + 0.005),  # Set y-axis limits
    expand = c(0, 0)  # No auto-expansion
  ) +
  labs(y = "Simpson", x = "") +
  theme_minimal() + 
  stat_summary(fun.data = mean_sd, geom = "errorbar", width= 0.2, colour = "black") +
  theme(
    axis.text.x = element_text(face = "bold", size = 17, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"), 
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  stat_pvalue_manual(
    data = sig_dunn,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )

dev.off()


############################ Test for normality
shapiro.test(df$Simpson)
# Result
#Shapiro-Wilk normality test

#data:  df$Simpson
#W = 0.69426, p-value < 2.2e-16
shapiro.test(df$Observed) 
# Result
#data:  df$Observed
#W = 0.91248, p-value = 1.29e-08
shapiro.test(df$Shannon)  
#Result
#data:  df$Shannon
#W = 0.86258, p-value = 2.086e-11
shapiro.test(df$Chao1) 
#Result
#W = 0.91246, p-value = 1.287e-08


