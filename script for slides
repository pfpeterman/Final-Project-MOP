```
#!/bin/bash
module load BLAST
for FILE in fastq_files/*; do
        BASE=$(basename "$FILE" .fastq.gz)
        seqtk seq -A "$FILE" > "${BASE}.fasta"
done
cat SRR*.fasta > concatenated.fasta
makeblastdb -in concatenated.fasta -dbtype nucl -out database_all
blastn -query reference.fasta -db database_all -out results_carp.txt -outfmt 6
```
