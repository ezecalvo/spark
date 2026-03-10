#!/bin/bash
# to run -> bash /pi/athma.pai-umw/home/jesse.lehman-umw/scripts/splicing_rates/bamrenamer /pi/athma.pai/... (dir_in)

# If hisat3n-mapped and converted read sorter used for nascent read identification:
# Loop through all files starting with "unique_converted"
for file in unique_converted*; do
  # Get the current file name without the prefix
  new_name="${file#unique_}"
  
  # Get the new file name with the addition of the naming convention
  new_name="${new_name%R1.bam}R1_15m_rep1.bam"
  # Rename the file
  mv "$file" "$new_name"
done


# If STAR-mapped and bam2bakR used for nascent read identification:
#Loop through all files starting with "unique_converted"
#for file in *R1Aligned.sortedByCoord.out.bam; do
#  # Get the new file name without the STAR suffixes and with the addition of the naming convention
#  new_name="${file%R1Aligned.sortedByCoord.out.bam}R1_15m_rep1.bam"
#  # Rename the file
#  mv "$file" "$new_name"
#done
