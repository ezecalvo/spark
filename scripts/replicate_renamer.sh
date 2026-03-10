#!/bin/bash
# to run: bash /pi/athma.pai-umw/home/jesse.lehman-umw/scripts/splicing_rates/replicate_renamer.sh "old prefix" "new prefix" "old rep" "new rep"

#old_prefix="converted_MDM7_0_R1"
#new_prefix="lps_0"
#old_rep="rep1"
#new_rep="rep1"

old_prefix=$1
new_prefix=$2
old_rep=$3
new_rep=$4

# === RENAME LOOP ===
for file in *"${old_prefix}"*"${old_rep}"*; do
    # Skip if no matching files
    [ -e "$file" ] || continue

    # Build new name
    new_file="${file//$old_prefix/$new_prefix}"
    new_file="${new_file//$old_rep/$new_rep}"

    # Rename the file
    echo "Renaming: $file → $new_file"
    mv "$file" "$new_file"
done

