####################################################################################
# Title: Heatmap (presence-absence) generation of Antibiotic Resistance Genes (ARGs)
# Author: [Ramani Shyam]
#
# Description:
# This script reads percent values of ARGs present different genus. Data are melted 
# to long format, and visualized as heatmaps using ggplot2. Each plot is exported as 
# a high-resolution PNG.
####################################################################################

# Loading libraries
library(ggplot2)
library(reshape2)
library(patchwork)
library(readxl)

sample_sizes <- c(8, 11, 29, 14, 69, 40, 33, 25, 31, 15, 22)  
y_labels <- paste(names(sample_sizes), "(", sample_sizes, ")", sep = "")

my_colors <- c("lightyellow", "orange", "violet", "black")

#Enzymatic ----
enzymatic <- read_excel("Enzymatic_ARGs_percent.xlsx")

melted_data_enzymatic <- reshape2::melt(enzymatic, id.vars = "ARGs",
                              variable.name = "Genus",
                              value.name = "Percent")

melted_data_enzymatic$ARGs <- factor(melted_data_enzymatic$ARGs,
                                     levels = unique(enzymatic$ARGs))

# Plotting heatmap using ggplot2
enzymatic_genus <- ggplot(melted_data_enzymatic, aes(x = ARGs, y = Genus, fill = Percent)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = my_colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16,face = "bold", color = "black"),
        axis.title.y = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, color = "black", face = "bold.italic"), legend.margin = margin(t = 50)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold", margin = margin(b = 15))) +
  labs(x = "ARGs (Enzymatic inactivation)", y = "Genus", fill = "% of genomes") +
  scale_y_discrete(labels = function(x) {paste(x, y_labels, sep = " ")})

# ARG and Class vectors
ARG_enzymatic <- c("aac(2')-Ia", "aac(3)", "aac(6')", "aac(6')-Ie-aph(2'')-Ia", "ant(2'')-Ia", "ant(3'')-IIa", "ant(6)-Ia", "ant(9)-Ia", "aph(3'')-Ib", "aph(3')", "aph(6)-Id", "aadA", "aad(6)", "blaACT-25", "blaADC", "amp_Ec", "blaAQU-3", "blaCARB-3", "blacphA", "blaCMH-1", "blaCMY", "blaCTX", "blaDHA", "blaDIM-1", "blaMIR", "blaMOX-7", "blaNDM", "blaOKP", "blaOXA", "blaOXY-1-1", "blaPDC", "blaSHV", "blaTEM", "blaVEB", "ereA", "erm", "mphA", "mphB", "mphE", "mphG", "tetX", "crpP", "ble", "lnuF", "sat-1", "sat-4", "cat", "fosA", "arr")
Class_enzymatic <- c("Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Beta-lactam", "Macrolide", "Macrolide", "Macrolide", "Macrolide", "Macrolide", "Macrolide", "Tetracycline", "Fluoroquinolone", "Glycopeptide", "Lincosamide", "Nucleoside", "Nucleoside", "Phenicol", "Phosphonic acid", "Rifamycin")
arg_groups_enzymatic <- data.frame(ARG = ARG_enzymatic, Drug_Class = Class_enzymatic)
arg_groups_enzymatic$ARG <- factor(arg_groups_enzymatic$ARG, levels = arg_groups_enzymatic$ARG)

# Annotation strip for ARG groups
annotation_strip_enzymatic <- ggplot(arg_groups_enzymatic, aes(x = ARG, y = 1, fill = Drug_Class)) +
  geom_tile(height = 0.1) +
  scale_fill_manual(values = c("Aminoglycoside" = "brown", "Beta-lactam" = "lightgreen", "Macrolide" = "lightblue", "Tetracycline" = "lightsalmon", "Fluoroquinolone" = "violet", "Nucleoside" = "yellow", "Rifamycin" = "darkgreen", "Phenicol" = "purple" , "Lincosamide" = "darkblue", "Phosphonic acid" = "pink", "Glycopeptide" = "orange")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.margin = margin(b = 10),
        legend.position = "top",
        theme(aspect.ratio = 0.05)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "Drug Class")

# Combining & saving the heatmap and annotation strip 
combined_plot_enzymatic <- annotation_strip_enzymatic + enzymatic_genus + plot_layout(heights = c(0.03, 1))
png(filename = "enzymatic_genus.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_enzymatic)
dev.off()

#Target ----
target <- read_excel("Target_ARGs_percent.xlsx")

melted_data_target <- reshape2::melt(target, id.vars = "ARGs",
                                        variable.name = "Genus",
                                        value.name = "Percent")

melted_data_target$ARGs <- factor(melted_data_target$ARGs,
                                     levels = unique(target$ARGs))

# Plotting heatmap using ggplot2
target_genus <- ggplot(melted_data_target, aes(x = ARGs, y = Genus, fill = Percent)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = my_colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16,face = "bold", color = "black"),
        axis.title.y = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size =12, angle = 0, hjust = 1, color = "black", face = "bold.italic"), legend.margin = margin(t = 50)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold", margin = margin(b = 15))) +
  labs(x = "ARGs (Target modification)", y = "Genus", fill = "% of genomes") +
  scale_y_discrete(labels = function(x) {paste(x, y_labels, sep = " ")})

# ARG and Class vectors
ARG_target <- c("armA", "rmt", "msrE", "msrC", "qnrA1", "qnrB", "qnrD1", "qnrS1", "qnrVC", "tetM", "dfr", "van", "lsaE", "mcr-9", "ugd", "arnA", "pmrF", "bacA", "basS", "eptA", "sul")
Class_target <- c("Aminoglycoside", "Aminoglycoside", "Macrolide", "Macrolide", "Fluoroquinolone", "Fluoroquinolone", "Fluoroquinolone", "Fluoroquinolone", "Fluoroquinolone", "Tetracycline", "Diaminopyrimidine", "Glycopeptide", "Lincosamide", "Peptide", "Peptide", "Peptide", "Peptide", "Peptide", "Peptide", "Peptide", "Sulfonamide")
arg_groups_target <- data.frame(ARG = ARG_target, Drug_Class = Class_target)
arg_groups_target$ARG <- factor(arg_groups_target$ARG, levels = arg_groups_target$ARG)

# Annotation strip for ARG groups
annotation_strip_target <- ggplot(arg_groups_target, aes(x = ARG, y = 1, fill = Drug_Class)) +
  geom_tile(height = 0.1) +
  scale_fill_manual(values = c("Aminoglycoside" = "maroon", "Macrolide" = "lightblue", "Fluoroquinolone" = "violet", "Diaminopyrimidine" = "yellow", "Glycopeptide" = "orange", "Lincosamide" = "purple", "Peptide" = "darkblue", "Sulfonamide" = "pink", "Tetracycline" = "lightsalmon")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.margin = margin(b = 10),
        legend.position = "top",
        theme(aspect.ratio = 0.05)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "Drug Class")

# Combining & saving the heatmap and annotation strip
combined_plot_target <- annotation_strip_target + target_genus + plot_layout(heights = c(0.03, 1))
png(filename = "target_genus.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_target)
dev.off()

#Efflux ----
efflux <- read_excel("Efflux_ARGs_percent.xlsx")

melted_data_efflux <- reshape2::melt(efflux, id.vars = "ARGs",
                                        variable.name = "Genus",
                                        value.name = "Percent")

melted_data_efflux$ARGs <- factor(melted_data_efflux$ARGs,
                                     levels = unique(efflux$ARGs))

# Plotting heatmap using ggplot2
efflux_genus <- ggplot(melted_data_efflux, aes(x = ARGs, y = Genus, fill = Percent)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = my_colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16,face = "bold", color = "black"),
        axis.title.y = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, color = "black", face = "bold.italic"), legend.margin = margin(t = 50)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold", margin = margin(b = 15))) +
  labs(x = "ARGS (Efflux)", y = "Genus", fill = "% of genomes") +
  scale_y_discrete(labels = function(x) {paste(x, y_labels, sep = " ")})

# ARG and Class vectors
ARG_efflux <- c("bae", "cpxA", "emrE_Pa", "kdpE", "ompK37_Kp", "abeS", "amvA_Ab", "emrE_Ec", "mefC", "ade", "mdfA_Ec", "tet(A)", "tet(B)", "tet(C)", "tet(D)", "tet(E)", "tet(G)", "tet(J)", "tet(L)", "tetU", "abaQ_Ab", "abeM", "qepA4", "efmA", "emr", "acr", "acrA_Enc", "acrA_Ec", "acrA_Kp", "armR", "axy", "crp", "cpxR_Pa", "evg", "gad", "hns", "kpn_Kp", "marA", "mdt", "mex", "mux", "opm", "opr", "oqx", "pmpM", "ramA", "soxR_Pa", "tolC", "abaF_Ab", "bcr-1", "cmIA", "floR", "msbA", "qacH", "triA", "triB", "triC", "yojI")
Class_efflux <- c("Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Aminoglycoside", "Beta-lactams", "Macrolide", "Macrolide", "Macrolide", "Macrolide", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Tetracycline", "Fluoroquinolone", "Fluoroquinolone", "Fluoroquinolone", "Macrolide-fluoroquinolone", "Tetracycline-fluoroquinolone", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Multidrug", "Others", "Others", "Others", "Others", "Others", "Others", "Others", "Others", "Others", "Others")
arg_groups_efflux <- data.frame(ARG = ARG_efflux, Drug_Class = Class_efflux)
arg_groups_efflux$ARG <- factor(arg_groups_efflux$ARG, levels = arg_groups_efflux$ARG)

# Annotation strip for ARG groups
annotation_strip_efflux <- ggplot(arg_groups_efflux, aes(x = ARG, y = 1, fill = Drug_Class)) +
  geom_tile(height = 0.1) +
  scale_fill_manual(values = c("Aminoglycoside" = "maroon", "Beta-lactam" = "lightgreen", "Macrolide" = "lightblue", "Tetracycline" = "lightsalmon", "Fluoroquinolone" = "violet", "Macrolide-fluoroquinolone" = "orange", "Tetracycline-fluoroquinolone" = "purple", "Multidrug" = "pink",  "Others" = "darkgreen")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.margin = margin(b = 10),
        legend.position = "top",
        theme(aspect.ratio = 0.05)) +
  theme(legend.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "Drug Class")

# Combine & saving the heatmap and annotation strip
combined_plot_efflux <- annotation_strip_efflux + efflux_genus + plot_layout(heights = c(0.03, 1))
png(filename = "efflux_genus.png", width = 20, height = 10, units = "in", res = 900)
print(combined_plot_efflux)
dev.off()