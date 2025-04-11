## Final-Project-MOP
final project repository for Mars, Olivia, Paige 

"Detecting Asian Carp Presence in the Great Lakes Region using Bioinformatics Tools"

our project question is how can bioinformatic tools be utilized to detect sturgeon presence within the Great Lakes?

process:
1. collect eDNA samples taken from the Great Lakes and nearby rivers such as the St. Lawerence River, The Grand River, and the Muskegon River.
2. use FASTQC to assess raw read quality of the files obtained.
3. Align eDNA reads to refernce sturgeon sequences using BLAST
4. Construct a map with QGIS.
5. Compare detection patterns across the Great Lakes Region. 

# Fork & clone repository before beginning

## Step 1: Create & Download Query

- Make a file for the accession numbers

```bash
vi srr_files.txt
```

- paste accession list:

```
SRR32172062
SRR32172063
SRR32172064
SRR32172065
SRR32172069
SRR32172073
SRR32172074
SRR32172075
SRR32172076
SRR32172077
SRR32543594
SRR32543595
SRR32543252
SRR32543253
SRR32543254
SRR32543256
SRR32543258
SRR32543259
SRR32543261
SRR32543260
SRR32543262
SRR32543263
SRR32543264
SRR32543265
SRR32543266
SRR32543267
SRR32543270
SRR32543269
SRR32543271
SRR32543272
SRR32543255
SRR32543592
SRR32543593
SRR32543592
SRR32543593
SRR32543555
SRR32543552
SRR32543553
SRR32543554
SRR32543564
SRR32172078
SRR32172079
SRR32172080
SRR32172081
SRR32172082
SRR32172083
SRR32172084
SRR32172085
SRR32172086
SRR32172087
SRR32172088
SRR32172089
SRR32172090
SRR32172091
SRR32172092
SRR32172093
SRR32172094
SRR32172095
SRR32172096
SRR32172103
SRR32172100
SRR32172097
SRR32172098
SRR32172099
SRR32172101
SRR32172102
SRR32172104
SRR32172110	
SRR32172105
SRR32172106
SRR32172113	
SRR32172107
SRR32172108
SRR32172109
SRR32172111
SRR32172112	
SRR32172114
SRR32543273
SRR32543274
SRR32543275
SRR32543348
SRR32543349
```

- Create script to download files from NCBI

```
vi download.sh
```

- paste the script into the file:

```
#!/bin/bash
module load sra-toolkit
while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR
   fastq-dump --gzip $SRR --split-files -O fastq_files/
done < srr_files.txt
```

- Execute the script

```
bash download.sh
```

## Step 2: Create Reference FASTA

- make reference file

```
vi reference.fasta
```

- paste the reference fasta (12S Gene)

```
>LC091587.1 Cyprinus carpio mitochondrial gene for 12S rRNA, partial sequence, specimen_voucher: UW:154625
CACCGCGGTTAGACGAGAGGCCCTAGTTGATATTACAACGGCGTAAAGGGTGGTTAAGGATAAACAAAAA
TAAAGTCAAATGGCCCCTTGGCCGTCATACGCTTCTAGGAGTCCGAAGCCCTAATACGAAAGTAACTTTA
ATAAACCCACCTGACCCCACGAAAGCTGAGAAA
```

## Step 3: BLAST query to reference

### Concatenating forward and reverse files into a single fasta file
- Create BLAST script 

```
vi concat.sh
```

- paste the script:

```
#!/bin/bash

mkdir -p concat

for fwd in *_1.fasta; do
    acc="${fwd%%_*}"
    rev="${acc}_2.fasta"

    if [[ -f "$rev" ]]; then
        cat "$fwd" "$rev" > "concat/${acc}.fasta"
        echo "Concatenated $acc"
    else
        echo "Missing reverse file for $acc"
    fi
done
```
### Running blast using each accession as a subject
- Create BLAST script 

```
vi blast.sh
```

- paste the script:

```
#!/bin/bash

mkdir -p blast_results

module load BLAST

# Loop over each subject file in the concat directory
for subject in concat/*.fasta; do
    # Extract accession from filename
    acc=$(basename "$subject" .fasta)

    # Make BLAST database (you could skip this if done once already)
    makeblastdb -in "$subject" -dbtype nucl -out "blastdb/${acc}"

    # Run BLAST
    blastn -query reference.fasta -db "blastdb/${acc}" \
           -out "blast_results/${acc}_vs_carp.txt" \
           -outfmt 6

  echo "Finished BLAST for $acc"

done
```
- NOTE: Default blast parameters above produce 500 hits. You can increase the number of hits by adding `-max-target-seq 1000` right after number 6 in the blastn command.
- NOTE: Identity threshold for fish species appears to be 99.45, lower identity may be questionable. Can you please run the blast script by adding this flag and see if you get many hits? -perc_identity 99.45

- execute script

```
bash blast.sh
```

combine results
```
cat blast_results/*.txt > results_carp.txt
```

## Step 4: Confirm Identity of Hits

- Extract hits from results

```
cut -f1 filtered_hits.txt | sort | uniq > matching_ids.txt
seqtk subseq concatenated.fasta carp_list.txt > carp_hits_seqs.fasta
```

- Open hits file

```
cat carp_hits_seqs.fasta
```
- select the hit with the lowest matching score
- Navigate to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
- Paste sequence and click BLAST
- Confirm the top hit is Cyprinus carpio

## Normalize hits
1. Cound number of reads per sample
```
vi read_count.sh
```
- Type I and paste:
```
#!/bin/bash
mkdir -p read_counts

for file in concat/*.fasta; do
    acc=$(basename "$file" .fasta)
    count=$(grep -c "^>" "$file")
    echo -e "${acc}\t${count}"
done > read_counts/read_counts.tsv
```
- Run read_count.sh script:
```
bash read_count.sh
```
Lets run the following in the terminal

```
cut -f1 blast_results/all_blast_results.tsv | sort | uniq -c | awk '{print $2"\t"$1}' > blast_hit_counts.tsv

awk 'FNR==NR{a[$1]=$2; next} ($1 in a){printf "%s\t%s\t%s\t%.2f\n", $1, a[$1], $2, (a[$1]/$2)*1000000}' \
blast_hit_counts.tsv read_counts/read_counts.tsv > carp_abundance_HPMR.tsv

sed -i '1iaccession\thit_count\tread_count\tHPMR' carp_abundance_HPMR.tsv

```


## Step 5: Create Figure

- Push results to repository
- Download .txt file 
- Make a new project in R Studio & set working directory to location of .txt file
- Paste and excecute the R script:

```
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
```

## Step 6: Create Map
