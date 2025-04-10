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

## Step 2: BLAST query to reference

- Create BLAST script 

```
vi blast.sh
```

- paste the script:

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

- execute script

```
bash blast.sh
```

## Step 3: Confirm Identity of Hits

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

## Step 4: Create Chart

- Push results to repository
- Download .txt file 
- Make a new project in R Studio & set working directory to location of .txt file
- Paste and excecute the R script:

```

```
