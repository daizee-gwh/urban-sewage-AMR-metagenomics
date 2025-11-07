####################################################################################
# Title: Heatmap (presence-absence) generation of Antibiotic Resistance Genes (ARGs)
# Author: [Ramani Shyam]
#
# Description:
# This script reads ARG presence-absence binary data from Excel sheets corresponding to
# different resistance mechanisms (Enzymatic inactivation, Target modification,
# and Efflux). Data are scaled, melted to long format, and visualized as
# heatmaps using ggplot2. Each plot is exported as a high-resolution PNG.
####################################################################################

# Loading libraries
library(ggplot2)
library(reshape2)
library(patchwork)
library(readxl)

#Enzymatic ----
Enzyme <- read_excel("Enzyme.xlsx")
Enzyme_scaled <- Enzyme
Enzyme_scaled[, c(3:307)] <- scale(Enzyme_scaled[, 3:307])

Enzyme$ARGs <- factor(x = Enzyme$ARGs,
                      levels = Enzyme_scaled$ARGs, 
                      ordered = TRUE)

enzymatic_305 <- ggplot(melt(Enzyme), aes(variable, ARGs, fill = Class, alpha = value)) + 
  geom_tile(colour = "gray50", linewidth = 0.1) +
  scale_alpha_identity(guide = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = NULL, hjust = NULL, color = "black", size = 8, face = "bold"), 
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"), plot.margin = margin(1, 10, 10, 10, "mm")) +
  theme(legend.text = element_text(size = 12, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  labs(x="Isolates", y="ARGs (Enzymatic inactivation)", fill = "Drug Class") 

#Target ----
Target <- read_excel("Target.xlsx")
Target_scaled <- Target
Target_scaled[, c(3:307)] <- scale(Target_scaled[, 3:307])

Target$ARGs <- factor(x = Target$ARGs,
                      levels = Target_scaled$ARGs, 
                      ordered = TRUE)

target_305 <- ggplot(melt(Target), aes(variable, ARGs, fill = Class, alpha = value)) + 
  geom_tile(colour = "gray50", linewidth = 0.1) +
  scale_alpha_identity(guide = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_text(angle = NULL, hjust = NULL, color = "black", size = 8, face = "bold"), 
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"), plot.margin = margin(1, 10, 10, 10, "mm")) +
  theme(legend.text = element_text(size = 12, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  labs(x="Isolates", y="ARGs (Target modification)", fill = "Drug Class") 

#Efflux ----
Efflux <- read_excel("Efflux.xlsx")
Efflux_scaled <- Efflux
Efflux_scaled[, c(3:307)] <- scale(Efflux_scaled[, 3:307])

Efflux$ARGs <- factor(x = Efflux$ARGs,
                      levels = Efflux$ARGs, 
                      ordered = TRUE)

efflux_305 <- ggplot(melt(Efflux), aes(variable, ARGs, fill = Class, alpha = value)) + 
  geom_tile(colour = "gray50", linewidth = 0.1) +
  scale_alpha_identity(guide = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_text(angle = NULL, hjust = NULL, color = "black", size = 8, face = "bold"), 
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"), plot.margin = margin(1, 10, 10, 10, "mm")) +
  theme(legend.text = element_text(size = 12, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  labs(x="Isolates", y="ARGs (Efflux)", fill = "Drug Class")

#annotation strip ----
Isolates <- c("C22S2B2", "C29S2B3", "C33S6B1", "C30S2B3_1", "UPS4B2", "C50S2B2", "C30S2B3_2", "C30S2B2", "CG1S1AB1", "C27S6B3", "C11S6B3", "C49S5B1", "C19S5B3", "C19S5B4", "C19S3B4", "C34S3B1", "C48S5B2", "KOL1AB1_3", "CG1S1B1", "KOL2B1", "C24S6B1", "C26S7B2", "C33S7B1", "C34S7B1", "JH1S1B5", "KOL1S3B1", "C22S7B1", "C28S6B1", "GUW4B1", "C24S4B1", "C24S6B2", "C24S5B1", "C26S7B3", "C31S3B1_1", "C36S3B1", "C36S7B1", "C34S4B1_2", "C46S2B1", "C43S6B1", "C48S3B2", "UPS6B3", "C31S3B1_2", "C45S2AB3", "UP1S5AB2", "DED_AB2", "HDW_AB1", "C29S4B1", "UPS3B2", "C35S4B1", "C36S3B2", "C36S6B1", "C41S6B2", "C47S1B2", "C41S1AB2", "C3S5B3", "C19S4B2", "C36S4B2", "C36S4B1", "GUWS5AB2", "C42S3B1", "C42S1B1", "RSK_1", "C7S5B2", "C8S5B2", "C48S5B4", "C49S4B1", "C12S5B2", "C13S2B1", "C3S5B2", "C6S3B2", "C13S3B1", "KOL1S1B2", "KOL1S2B2", "CG1S1B2", "CG1S1B3", "C19S3B3", "C19S1B3", "C19S6B1", "C20S1B2", "C20S4B1", "C21S6B2", "C22S1B2", "C22S2B1", "C22S4B1", "C22S6B1", "C23S1B2", "C23S2B3", "C23S3B2", "C23S4B3", "C23S6B2", "C23S7B1", "C24S3B2", "C25S2B2", "C25S3B2", "C25S4B3", "C26S1B1", "C26S5B2", "C26S6B2", "C28S2B1", "C28S3B2", "C28S4B3", "C28S5B2", "C28S6B3", "C29S5B1", "C30S5B1", "C30S6B1", "C33S2B1", "C33S6B2", "C35S1B1", "C35S2B1", "C35S2B2", "C35S3B2", "C35S4B2", "C30S1B2", "C31S2B1", "C32S6B2", "C34S7B2", "C34S3B2", "C36S6B2", "C35S6B2", "C35S7B2", "C44S7B2", "C49S1B1", "UPS2B4", "UPS1B1", "UPS5B1", "UPS10B1", "DED-1", "HDW-2", "JH1S2B2", "C49S1B2", "C43S5B1_1", "CG3S4AB1", "C41S6AB1", "C43S5B1_2", "DED_AB1", "C49S6AB1", "C49S6B2", "C19S6B2", "C20S2B2", "C20S3B2", "C21S1B1", "C22S3B1", "C23S2B2", "C24S2B1", "C25S4B2", "C25S5B2", "C25S6B2", "C26S2B1", "C27S5B3", "C29S4B2", "C33S2B2", "C35S5B1", "C31S4B1", "C31S5B1", "C34S6B2", "C34S5B1", "C32S4B2", "C41S2B2", "C37S3B1", "C43S1B1", "C45S6B1", "C43S3B1", "C45S1B1", "C48S4B2", "UPS4B4", "UPS7B1", "UPS9B1", "UPS10B2", "C49S4B2", "C35S2AB1", "C19S3B1", "C26S3B2", "C26S5B3", "C28S3B1", "C40S4B2", "C50S4B3", "C2S3B1", "C4S2B1", "C7S4B1", "C10S4B1", "C10S6B1", "C11S1B1", "C11S5B1", "C11S6B1", "C11S7B1", "C10S5AB1", "C10S7AB1", "C11S2AB1", "C12S2AB1", "C12S4AB1_1", "C15S9B3", "C17S5B1", "C18S3B1", "HR1S1B1", "C19S3B2", "C19S4B1", "C19S5B2", "C12S4AB1_2", "C14S6AB1", "C23S4B2", "C49S2B3_1", "C44S3B1", "C48S1B3", "UPS2B2", "UPS4B3", "CG3S7B2", "C49S2B3_2", "C50S1B2", "UP1S6AB2", "C1S1B2", "C1S3B1", "C1S4B1", "C11S3B2", "C11S6AB2", "C20S3B1", "C21S5B1", "C15S2AB1", "C15S9AB1", "C33S1B1", "C42S4B1", "C42S5B1", "C47S6B3", "UPS3B3", "C37S1AB1", "C42S2AB1_1", "UP1S3AB1", "UP1S4AB1", "C45S6AB2", "C29S6B2", "C42S2AB1_2", "C45S3AB3", "UPS8B1", "C24S7B1", "RSK_AB1", "C5S1B1", "C6S4B1", "C7S6B1", "C8S5B1", "C9S2B1", "C9S5B2", "C11S4B2", "C11S5B2", "C11S7B2", "C12S2B2", "C12S3B1", "C15S3B3", "C16S5B2", "C17S3B2", "C10S4B2", "C2S2B1", "C5S3B2", "C9S5B1", "C12S2B1", "C12S5B1", "C13S4B1", "C23S3B1", "C34S4B1_1", "C48S3B1", "C50S3B2", "C4S1B1", "C11S6B2", "C1S4B2", "C4S4B1", "C15S9B2", "C17S6B2", "C20S2B1", "C23S5B1", "C24S5B2", "C41S2B1", "C39S6B2", "C43S5B3", "KOL1S1B1", "KOL1S3B2", "KOL1S3B3", "C22S6B2", "C29S7B1", "JH1S10B1", "KOL1S2B1", "C31S6B2", "CG3S1B1", "C10S1AB2", "C10S2AB2", "C10S6AB2", "C10S7AB2", "C11S3AB1", "C11S4AB1", "C11S5AB1", "C12S1AB1", "C12S3AB1", "C12S4AB2", "C14S1AB1", "C14S7AB2", "C16S1AB1", "C25S2AB2", "C25S4AB2", "C26S2AB1", "C26S5AB1", "C26S7AB1", "CG2S1AB1", "JHS3AB2", "CG3S6AB1", "C25S1AB1", "C38S1B2", "C38S3B2")

Species <- c(rep("A. baumannii", 8),
             rep("A. caviae", 2), "A. dhakensis", "A. simiae", rep("A. veronii", 6), "Aeromonas sp.", 
             "Achromobacter sp.", "Alcaligenes sp.", "Alcaligenes sp.", "Alcaligenes sp.", "B. trematum", "C. aquatica",
             rep("C. freundii", 4), "C. portucalensis", rep("C. werkmanii", 18), rep("Citrobacter sp.", 6), 
             rep("E. cloacae", 5), rep("E. hormaechei", 3), rep("E. roggenkampii", 5), "Enterobacter sp.", 
             rep("E. coli", 69), "K. michiganensis", rep("K. pneumoniae", 33), rep("K. quasipneumoniae", 6),
             rep("M. morganii", 33), rep("P. mirabilis", 22), rep("P. penneri", 1), rep("P. vulgaris", 2), 
             rep("P. manganoxydans", 15), rep("P. rettgeri", 10), rep("P. stuartii", 2), rep("Providencia sp.", 4), 
             rep("P. aeruginosa", 6), rep("P. putida", 6), rep("Pseudomonas sp.", 3),
             rep("E. faecium", 21), rep("E. raffinosus", 1), rep("L. fusiformis", 2))
samples <- data.frame(
  Isolates = factor(Isolates, levels = Isolates, ordered = TRUE),
  Species = factor(Species)
)
# Creating the annotation strip for the isolates
annotation_strip <- ggplot(samples, aes(x = Isolates, y = 1, fill = Species)) +
  geom_tile(height = 0.05) +
  scale_fill_manual(values = c("A. baumannii" = "#FF0000", "E. faecium" = "#E0115F", "E. coli" = "#B7410E", "A. caviae" = "#E30B5D", "A. simiae" = "#272124", "A. dhakensis" = "#4169E1", "A. veronii" = "#00755E", "Aeromonas sp." = "#4B0082", "Alcaligenes sp." = "#826644", "C. werkmanii" = "#00CCCC", "C. freundii" = "#AA98A9", "Achromobacter sp." = "#B76E79", "C. portucalensis" = "#7851A9", "E. cloacae" = "#CD5909", "B. trematum" = "#FF355E", "E. hormaechei" = "#004225", "E. roggenkampii" = "#A45A52", "Enterobacter sp." = "#838996", "C. aquatica" = "#BC8F8F", "K. pneumoniae" = "#8D4E36", "K. quasipneumoniae" = "#E3256B", "K. michiganensis" = "#522D80", "P. aeruginosa" = "#FBAB60", "P. putida" = "#6D3B24", "Pseudomonas sp." = "#8A7F80", "M. morganii" = "#9B111E", "P. manganoxydans" = "#BB6528", "P. rettgeri" = "#480607", "Providencia sp." = "#32174D", "P. mirabilis" = "#7F7F7F", "Citrobacter sp." = "#D27D46", "P. vulgaris" = "#8DA8CC", "P. penneri" = "#00FFC0", "E. raffinosus" = "#4B86B4", "L. fusiformis" = "#6B4226", "P. stuartii" = "#9B1C31")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.margin = margin(t = 10),
        legend.position = "top",
        theme(aspect.ratio = 0.05)) +
  theme(legend.text = element_text(size = 10, color = "black", face = "bold.italic")) +
  theme(legend.title = element_text(size = 12, color = "black", face = "bold")) +
  guides(fill = guide_legend(ncol = 9, bycol = TRUE))

# Combining the heatmap and annotation strip
combined_plot_enzymatic <- annotation_strip + enzymatic_305 + plot_layout(heights = c(0.03, 1)) 
combined_plot_target <- annotation_strip + target_305 + plot_layout(heights = c(0.03, 1))
combined_plot_efflux <- annotation_strip + efflux_305 + plot_layout(heights = c(0.03, 1))

# Saving the images
png(filename = "enzymatic_arg.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_enzymatic)
dev.off()

png(filename = "target_arg.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_target)
dev.off()

png(filename = "efflux_arg.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_efflux)
dev.off()