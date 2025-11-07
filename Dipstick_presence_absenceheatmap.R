#Author:Daizee Talukdar
# Load necessary library
library(ggplot2)
library(reshape2)
getwd()
# Input data
# Remove the first column (Dipstick_PCR_No) as it's an identifier
data <- as.data.frame(read.csv("Table_S2_R1.csv"))
data_matrix <- data[, -1]
rownames(data_matrix) <- data$Samples
# Reshape the data to long format
data_long <- melt(data, varnames = c("Sample", "Gene"), value.name = "Presence")
png(filename = "image3.png", 
    width = 10, height = 14, units = "in", 
    bg = "white",  res = 600)
# Create the heatmap
ggplot(data_long, aes(x = variable, y = Samples, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "gray", "1" = "darkblue"), name = "Presence") +
  labs(
    #title = "Heatmap of Binary Values",
    x = "Genes",
    y = "Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
dev.off()
