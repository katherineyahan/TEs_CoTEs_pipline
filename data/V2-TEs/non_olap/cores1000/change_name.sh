#!/bin/bash

# Directory path
directory="/cluster/home/t124771uhn/data/V2-TEs/non_olap/cores1000"

# Loop through files and rename
for file in "$directory"/*.bed; do
    new_filename="${file/_cores/}"
    mv "$file" "$new_filename"
    echo "Renamed $file to $new_filename"
done

echo "Renaming complete."

