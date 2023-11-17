#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yahan.zhang@mail.utoronto.ca
#SBATCH --job-name=generate_table
#SBATCH --output=TEs_CoTEs_pipline_1117.o
#SBATCH --error=TEs_CoTEs_pipline_1117.e

module load python3/3.10.9


# Assuming your folders are named cores1000 and TEs
folder1="/cluster/projects/lupiengroup/People/yahan/rotation_project/TEs_CoTE_pipline/data/V2-TEs/non_olap/TEs"
folder2="/cluster/projects/lupiengroup/People/yahan/rotation_project/TEs_CoTE_pipline/data/V2-TEs/non_olap/cores1000"

# Iterate through files in the first folder
for file1 in "$folder1"/*; do
    # Extract the file name without the path and extension
    file1_name=$(basename -- "$file1")
    file1_base="${file1_name%.*}"

    # Check if the corresponding file exists in the second folder
    file2="$folder2/$file1_base.bed"
    if [ -e "$file2" ]; then
        # Perform actions on the two files with the same name
        echo "Processing files: $file1 and $file2"

        # Add your actions here
	python TEs_CoTEs_script_1117.py $file1 $file2
    else
        echo "File $file2 does not exist in $folder2"
    fi
done
