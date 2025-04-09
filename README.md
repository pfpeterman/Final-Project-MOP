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

## Step 1: Create & Download Query

# 1. Make a file for the accession numbers

```bash
vi srr_files.txt
```

paste accession list, one number per line

# 2. Create script to download files from NCBI

```
vi download.sh
```

paste the script into the file:

```
#!/bin/bash
module load sra-toolkit
while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR
   fastq-dump --gzip $SRR --split-files -O fastq_files/
done < srr_files.txt
```

# 3. Execute the script

```
bash download.sh
```

## Step 2: BLAST query to reference

# 1. create BLAST script 

```
vi blast.sh
```

paste the script:

```
#!/bin/bash

module load BLAST

#converting fastq to fasta
for FILE in fastq_files/*; do

        BASE=$(basename "$FILE" .fastq.gz)
        seqtk seq -A "$FILE" > "${BASE}.fasta"
done

cat SRR*.fasta > concatenated.fasta

#making subject database (run only once)
makeblastdb -in concatenated.fasta -dbtype nucl -out database_all

#run BLAST
blastn -query reference.fasta -db database_all -out results_carp.txt -outfmt 6
```

## Step 3: Collect Data

# 1. Extract hits from results

```
cut -f1 filtered_hits.txt | sort | uniq > matching_ids.txt
seqtk subseq concatenated.fasta carp_list.txt > carp_hits_seqs.fasta
```

# 2. Conver blast results into an abundance table
