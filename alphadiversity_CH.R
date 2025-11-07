#Author:Daizee Talukdar

# load the library
library(ggplot2)
library(ggpubr)
library(dplyr)

library(writexl)

# Write your dataframe (df) to an Excel fil
df<-read.csv("Diversity_metadatasamples_updated.csv", header = TRUE)

p3 <- ggplot(df, aes(x = Months, y = Observed, fill = settings)) +
  geom_boxplot(width = 0.5, outlier.size = 2, outlier.shape = 21) +  # Boxplot with adjusted width and outlier shape
  
  theme_classic() +
  labs(x = NULL, y = "Observed") + theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey80")
  ) + theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
     
df$Observed <- as.numeric(df$Observed)

stat.test2 <- df %>%
  group_by(Months) %>%
  filter(n_distinct(settings) == 2) %>%     # keep only months with both groups
  wilcox_test(Observed ~ settings) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

stat.test2



p3 <- p3 + stat_pvalue_manual(stat.test2, hide.ns = TRUE)

ggsave("Observed.png", plot = p3, width = 10, height = 5, dpi = 600)







