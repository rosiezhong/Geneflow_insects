#!/bin/bash
module load Anaconda3/2022.05
source activate ~/genomes

# Exit if there is any error
set -e

# Get species name from the first argument
species="$1"

# Ensure species name is provided
if [ -z "$species" ]; then
    echo "Usage: $0 <species_name>"
    exit 1
fi

seqkit grep -n -r -p "chromosome" "${species}2.1.fasta.gz" > "${species}2.1.chromosome.fasta"
seqkit grep -n -r -p "contig|organelle" "${species}2.1.fasta.gz" > "${species}2.1.contig.fasta"
seqkit fx2tab --length --name --header-line "${species}2.1.chromosome.fasta" > "${species}2.1.sum.txt"
sum=$(awk '{sum+=$8} END {print sum}' "${species}2.1.sum.txt")
y=$(zcat "${species}2.2.fasta.gz" | grep ">" | wc -l)
z=$(cat "${species}2.1.contig.fasta" | grep ">" | wc -l)
window_size=$(echo "($sum / ($y - $z))" | bc -l)
window_size=$(printf "%.0f" "$window_size")
seqkit sliding -s "$window_size" -W "$window_size" "${species}2.1.chromosome.fasta" -o "${species}2.1.chromosome.split.fasta"
cat "${species}2.1.chromosome.split.fasta" "${species}2.1.contig.fasta" > "${species}2.1.split.fasta"
gzip "${species}2.1.split.fasta"
rm "${species}2.1.fasta.gz"
mv "${species}2.1.split.fasta.gz" "${species}2.1.fasta.gz"
rm "${species}2.1.chromosome.fasta" "${species}2.1.contig.fasta" "${species}2.1.sum.txt" "${species}2.1.chromosome.split.fasta"
unset sum
unset y
unset z
unset window_size

# Print success message
echo "Process completed successfully for $species"