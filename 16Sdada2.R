# To generate amplicon sequence variants from raw sequencing data. (Paired end reads)

# Load the following library
library(dada2); packageVersion("dada2")
path <- "/file"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#plot qulatiy
png(filename = "F.png",width = 12,height = 10,units = "in",res = 600)
plotQualityProfile(fnFs[1:2])
dev.off()
png(filename = "R.png",width = 12,height = 10,units = "in",res = 600)
plotQualityProfile(fnRs[1:2])
dev.off()
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#learning errrors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
png(filename = "Error_F.png",width = 12,height = 10,units = "in",res = 600)
plotErrors(errF, nominalQ=TRUE)
dev.off()
png(filename = "Error_R.png",width = 12,height = 10,units = "in",res = 600)
plotErrors(errR, nominalQ=TRUE)
dev.off()
#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#merging paried end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#constructing sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#write(seqtab, "seqtab_RUN84.txt")
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#removal of chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#save file to check sequences
write.csv(seqtab.nochim, "seq.csv")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# Track reads through pipeline
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
# Save the file to reported as results 
write.csv(track,  "read.csv")
# ASSIGN TAXONOMY
taxa <- assignTaxonomy(seqtab.nochim, "/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1

# preparing the ASV table (merging the taxonomic table and abundance table)
df1 <- t(otu_table(ps1))
df2 <- tax_table(ps1)
df3 <- cbind(df2, df1)
# Save the file
write.csv(df3, "ASV.csv")

# Retain bacterial taxa
ps_bacteria <- subset_taxa(ps, Kingdom == "Bacteria")
#Filtering
physeq_fil = filter_taxa(ps_bacteria, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# Ampvis code
install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2", Ncpus = 96)

library(ampvis2)

# prepare ampvis2 phyloseq
datad <- amp_load(
  otutable = "asv.csv",
  metadata = "metadata.csv",
  taxonomy = "taxonomy.csv"
)
#heatmap for Fig 3
amp_heatmap(
  datad,
  group_by = "sample",
  facet_by = "Group",
  tax_aggregate = "Genus",
  tax_add = "Species",
  tax_show = 5,
  color_vector = c("royalblue4","whitesmoke", "darkred"),
  plot_colorscale = "sqrt",
  plot_values = FALSE) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")
# heatmap for Fig 5
#final code for genus panindia new taxa names

amp_heatmap(
  data = datad_bitsp,
  group_by = "months",
  tax_show = 20,plot_colorscale = "log10",
  color_vector = brewer.pal (5, "Spectral"),
  tax_aggregate = "Genus",tax_add = "Phylum", min_abundance = 0.1,showRemainingTaxa = FALSE,
  plot_values = TRUE, plot_values_size = 2) + theme(axis.text.x = element_text(angle = 45, size=5, vjust = 1),
                                                    axis.text.y = element_text(size=9, face = "italic"),
                                                    legend.position="right")
# Calculate alpha diversity
# Assume your phyloseq object is called ps
alpha_div <- estimate_richness(physeq_fil)

