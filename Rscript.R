# Load required library
library(ggplot2)

# Load BLAST output from results_carp.txt
blast_data <- read.table("results_carp.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Assign column names based on BLAST outfmt 6
colnames(blast_data) <- c(
  "sseqid", "qseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

# Extract SRA accession (prefix before the dot in query ID)
blast_data$accession <- sub("\\..*", "", blast_data$qseqid)

# Count hits per accession
hit_counts <- as.data.frame(table(blast_data$accession))
colnames(hit_counts) <- c("Accession", "Hit_Count")

# Sort by hit count descending
hit_counts <- hit_counts[order(hit_counts$Hit_Count, decreasing = TRUE), ]

# Plot horizontal bar chart
ggplot(hit_counts, aes(x = reorder(Accession, Hit_Count), y = Hit_Count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Carp BLAST Hits per SRA Accession",
    x = "SRA Accession",
    y = "Number of Hits"
  ) +
  theme_minimal(base_size = 14)
