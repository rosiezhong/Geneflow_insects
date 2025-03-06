#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=150g
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH -o /Path/%x.o

VAR_DIR=
TREE_DIR=
species=
ratio=
end=
vcf_num=
SINGER_DIR=
WOR_DIR=
mkdir "$VAR_DIR/${species}.arg"
mkdir "$TREE_DIR/${species}.tskit"
touch "$VAR_DIR/${species}.singer_log.txt"
find "$VAR_DIR/${species}.phased.vcfs/" -type f | shuf -n "$vcf_num" | sed 's/\.vcf$//' > "$VAR_DIR/${species}_100.txt"

for line in $(cat "$VAR_DIR/${species}_100.txt"); do
    filename=$(basename "$line")
    echo "singer_master started for $filename"
    "$SINGER_DIR/singer-0.1.7/singer_master" -m 2.9e-9 -ratio "$ratio" -vcf "$line" -output "$VAR_DIR/${species}.arg/${filename}.arg" -start 0 -end "$end" -n 100 -thin 20
    if [ $? -ne 0 ]; then
        echo "Error with singer_master on $filename, skipping." >> "$VAR_DIR/${species}.singer_log.txt"
        rm -f -- "$VAR_DIR/${species}.arg"/*
        continue
    fi
    "$SINGER_DIR/singer-0.1.7/convert_to_tskit" -input "$VAR_DIR/${species}.arg/${filename}.arg" -output "$TREE_DIR/${species}.tskit/${filename}.tskit" -start 50 -end 99 -step 1
    if [ $? -ne 0 ]; then
        echo "Error with convert_to_tskit on $filename, skipping." >> "$VAR_DIR/${species}.singer_log.txt"
        rm -f -- "$VAR_DIR/${species}.arg"/*
        continue
    fi
    echo "Success for $filename." >> "$VAR_DIR/${species}.singer_log.txt"
    rm -f -- "$VAR_DIR/${species}.arg"/*
done
echo "Singer_master completed for $species"

python "$WOR_DIR/tree_type.py" "$TREE_DIR/${species}.tskit" 2 "$TREE_DIR/${species}.tree.txt"
if [[ $? -ne 0 ]]; then
    echo "Error: Generating tree types failed for $species"
    exit 1
fi
echo "Generating tree types completed for $species"
