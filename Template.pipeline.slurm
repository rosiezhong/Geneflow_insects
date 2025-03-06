#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=100g
#SBATCH --ntasks-per-node=16
#SBATCH --time=240:00:00
#SBATCH -o /data/zool-barralab/biol0249/working/output/%x.o

module load minimap2/2.24-GCCcore-11.3.0
module load HTSlib/1.14-GCC-11.2.0
module load SAMtools/1.16.1-GCC-11.3.0
module load R/4.4.0-gfbf-2023a

species_dir=
ALN_DIR=
VAR_DIR=
SCR_DIR=
WOR_DIR=

species=$(basename "$species_dir")
set -e  
set -o pipefail
chmod +x "$WOR_DIR/chr1.sh"
source "$WOR_DIR/chr1.sh"
ref="${species_dir}/${species}1.1.fasta.gz"

if [ ! -f "$ref" ]; then
    echo "Reference file $ref does not exist. Exiting."
    exit 1
fi

set -e  
set -o pipefail
chmod +x "$WOR_DIR/split.sh"
source "$WOR_DIR/split.sh"

for READS in 1.1 1.2; do
    minimap2 -ax asm5 --sam-hit-only --secondary=no -Y -R "@RG\tID:${species}${READS}\tSM:${species}${READS}" "$ref" "${species_dir}/${species}${READS}.fasta.gz" | samtools sort -@ 8 -m 4G -T "$ALN_DIR/temp_${species}${READS}" - > "$ALN_DIR/${species}${READS}.bam" && samtools index "$ALN_DIR/${species}${READS}.bam"
    if [ $? -ne 0 ]; then
    echo "Error doing alignment for ${species}${READS}. Exiting."
    exit 1
    fi
done 

for READS in 2.1 2.2; do
    minimap2 -ax asm20 --sam-hit-only --secondary=no -Y -R "@RG\tID:${species}${READS}\tSM:${species}${READS}" "$ref" "${species_dir}/${species}${READS}.fasta.gz" | samtools sort -@ 8 -m 4G -T "$ALN_DIR/temp_${species}${READS}" - > "$ALN_DIR/${species}${READS}.bam" && samtools index "$ALN_DIR/${species}${READS}.bam"
    if [ $? -ne 0 ]; then
    echo "Error doing alignment for ${species}${READS}. Exiting."
    exit 1
    fi
done 

echo "Alignment by minimap2 completed."

samtools merge -fr -o "$ALN_DIR/${species}1.merged.bam" "$ALN_DIR/${species}1.1.bam" "$ALN_DIR/${species}1.2.bam"
samtools merge -fr -o "$ALN_DIR/${species}2.merged.bam" "$ALN_DIR/${species}2.1.bam" "$ALN_DIR/${species}2.2.bam" && \
echo "Both bam files successfully merged."
callvariants2.sh in="$ALN_DIR/${species}1.merged.bam,$ALN_DIR/${species}2.merged.bam" ref="$ref" vcf="$VAR_DIR/${species}.vcf" ploidy=2 minreads=1 supplimentary=t

if [ $? -ne 0 ]; then
    echo "Error calling variants for $species. Exiting."
    exit 1
fi
echo "Variant calling by callvariants.sh completed."

rm "$SCR_DIR/individual_${species}1.merged.vcf.gz" "$SCR_DIR/individual_${species}2.merged.vcf.gz" "$ALN_DIR/${species}1.1.bam" "$ALN_DIR/${species}1.1.bam.bai" "$ALN_DIR/${species}1.2.bam" "$ALN_DIR/${species}1.2.bam.bai" "$ALN_DIR/${species}2.1.bam" "$ALN_DIR/${species}2.1.bam.bai" "$ALN_DIR/${species}2.2.bam" "$ALN_DIR/${species}2.2.bam.bai"

samtools depth "$ALN_DIR/${species}1.merged.bam" -o "$ALN_DIR/${species}1.depth.txt"
samtools depth "$ALN_DIR/${species}2.merged.bam" -o "$ALN_DIR/${species}2.depth.txt" 
echo "Bam files merged."

mkdir "$VAR_DIR/${species}.bedfiles"
mkdir "$VAR_DIR/${species}.vcfs"
mkdir "$VAR_DIR/${species}.phased.vcfs"
Rscript "$WOR_DIR/make_bedfiles.R" "$ALN_DIR/${species}1.depth.txt" "$ALN_DIR/${species}2.depth.txt" "$VAR_DIR/${species}.vcf" "$VAR_DIR/${species}.bedfiles/${species}"

if [ $? -ne 0 ]; then
    echo "Error in creation of bedfiles for $species. Exiting."
    exit 1
fi
echo "Creation of bedfiles for $species completed."

rm "$ALN_DIR/${species}1.depth.txt" "$ALN_DIR/${species}2.depth.txt"
bgzip "$VAR_DIR/${species}.vcf"
module load BCFtools/1.14-GCC-11.2.0
count=$(find "$VAR_DIR/${species}.bedfiles" -type f | wc -l)
tabix "$VAR_DIR/${species}.vcf.gz"

for INDEX in $(seq 1 "$count"); do
    bcftools view -Ov -R "$VAR_DIR/${species}.bedfiles/${species}${INDEX}.bed" "$VAR_DIR/${species}.vcf.gz" > "$VAR_DIR/${species}.vcfs/${species}${INDEX}.vcf"
    if [ $? -ne 0 ]; then
        echo "Error creating individual vcfs per region for $species. Exiting."
        exit 1
    fi
done
echo "Creation of individual vcfs per region for $species completed."

for INDEX in $(seq 1 "$count"); do
    sed -E 's/([01])\/([01])/\1|\2/g'  "$VAR_DIR/${species}.vcfs/${species}${INDEX}.vcf" > "$VAR_DIR/${species}.phased.vcfs/${species}${INDEX}.phased.vcf"
    if [ $? -ne 0 ]; then
        echo "Error phasing vcfs for $species. Exiting."
        exit 1
    fi
done
echo "vcfs for $species phased."

for INDEX in $(seq 1 "$count"); do
    variant_count=$(bcftools view "$VAR_DIR/${species}.phased.vcfs/${species}${INDEX}.phased.vcf" | grep -v '^#' | wc -l) 
    if [ "$variant_count" -le 100 ]; then
        rm "$VAR_DIR/${species}.phased.vcfs/${species}${INDEX}.phased.vcf" 
    fi
    if [ $? -ne 0 ]; then
        echo "Error removing vcfs with =<100 variants for $species. Exiting."
        exit 1
    fi 
done
echo "vcfs with =<100 variants for $species removed."
unset variant_count
unset count

if [ -z "$(find "$VAR_DIR/${species}.phased.vcfs" -maxdepth 1 -type f)" ]; then
    echo "Error: The directory ${species}.phased.vcfs is empty. Exiting."
    exit 1
else
    echo "Directory ${species}.phased.vcfs contains files. Proceeding."
fi

rm -rf "$VAR_DIR/${species}.bedfiles"
rm -rf "$VAR_DIR/${species}.vcfs"

echo "Proceeding to bgzip bam files."
bgzip "$ALN_DIR/${species}1.merged.bam"
bgzip "$ALN_DIR/${species}2.merged.bam"
# rm "$ALN_DIR/${species}1.merged.bam" "$ALN_DIR/${species}2.merged.bam" "$VAR_DIR/${species}.vcf.gz"
