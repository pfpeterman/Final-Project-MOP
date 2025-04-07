# Final-Project-MOP
final project repository for Mars, Olivia, Paige 

"Detecting Sturgeon Presence in the Great Lake Region using Bioinformatic Tools"

our project question is how can bioinformatic tools be utilized to detect sturgeon presence within the Great Lakes?

process:
1. collect eDNA samples taken from the Great Lakes and nearby rivers such as the St. Lawerence River, The Grand River, and the Muskegon River.
2. use FASTQC to assess raw read quality of the files obtained.
3. use Trimmomatic to trim out low quality data
4. Align eDNA reads to refernce sturgeon sequences using BLAST
5. Alternatively, use Kraken2 for taxonomic classification of metagenomic reads.
6. Extract matching sequences and align them using MAFFT.
7. Construct a phylogeneic tree using FASTtree.
8. Construct a map with QGIS.
9. Compare detection patterns across the Great Lakes Region. 


## Extracting hit labels from blast results

cut -f1 filtered_hits.txt | sort | uniq > matching_ids.txt
seqtk subseq concatenated.fasta carp_list.txt > carp_hits_seqs.fasta

## Conver blast results into an abundance table
