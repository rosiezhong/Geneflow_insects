#!/bin/bash
seqkit head -n 1 "${species_dir}/${species}1.1.fasta.gz" > "${species_dir}/${species}1.1chr1.fasta"
rm "${species_dir}/${species}1.1.fasta.gz"
bgzip "${species_dir}/${species}1.1chr1.fasta"
mv "${species_dir}/${species}1.1chr1.fasta.gz" "${species_dir}/${species}1.1.fasta.gz"

# Print success message
echo "Chr1 extracted successfully for $species"