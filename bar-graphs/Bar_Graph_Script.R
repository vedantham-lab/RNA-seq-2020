# Load libraries
library(tibble)
library(dplyr)
library(readxl)

# Load raw normalized count data
norm_counts <- read.csv("./Data/normalized_counts.csv",header = TRUE)
# List of genes to plot
genes_to_plot <- c("Isl1","Hcn4","Tbx3","Shox2","Bmp4","Tbx18",     #SAN>RA
                   "Tbx5","Mef2c","Gata4",                          #SAN~RA
                   "Nkx2-5","Nppa","Nppb","Bmp10","Gja1","Gja5")    #SAN<RA

# Sanity check
mean(genes_to_plot %in% norm_counts$Gene)
# [1] 1

# Get the data for the genes to plot
data_to_plot <- norm_counts[which(norm_counts$Gene %in% genes_to_plot),]
# Save the data and format in Excel
write.csv(data_to_plot, file="./Data/norm_count_data_for_bar_graph.csv", row.names = FALSE)

## Edit the data in Excel and format it for bar graph ##

# Read the formatted data back into R for prepararation of bar graph
data_to_plot <- read_excel("./Data/norm_count_data_for_bar_graph.xls",sheet = 2)

data_to_plot$Gene <- factor(data_to_plot$Gene, levels=unique(data_to_plot$Gene))
data_to_plot$Replicate <- factor(data_to_plot$Replicate, levels=unique(data_to_plot$Replicate))

# Plotting
library(ggplot2)
bar <- ggplot(data_to_plot, aes(x=Replicate, y=log2FC,fill=Status)) +
  geom_bar(stat='identity') +
  facet_wrap(~Gene,ncol = length(unique(data_to_plot$Gene))) +
  ylim(-11.5, 11.5) +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_hline(yintercept=-1, linetype="dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.text.x = element_text(face="bold.italic")) +
  theme(axis.title.x = element_blank()) +
  ylab("log(2) Normalized Counts PC/RACM") +
  scale_fill_discrete(breaks=c("PC>RACM", "PC<RACM", "PC~RACM"))
bar
ggsave(filename = "barplot_log2_norm_counts_PC_RACM.png", plot = bar, width = 30, height = 15, dpi = 300, units = "cm")
  
