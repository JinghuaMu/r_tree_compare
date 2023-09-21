#!/bin/bash

# Set the paths and variables
iqtree_loc="/home/timmu/miniconda3/envs/phylo/bin/iqtree2"
geneset_loc="../datasets/different_method/Squamates-Final_Alns/MAFFT-Auto-gt/"
outgroup_file="../datasets/different_method/Squamates-Final_Alns/Squamates-Taxon-Outgroup.txt"
species_tree_file="/home/timmu/project/mixmodel/datasets/different_method/phylogenic_tree/Squamates/Rooted_MAFFT_Auto_GapThreshold_Species_Tree_Scored.tre"

# Read the gene sets
gene_sets=( $(ls $geneset_loc) )
gene_sets=($(printf "%s\n" "${gene_sets[@]}" | sort))

echo "Current gene file location: $geneset_loc"
echo "Number of genes: ${#gene_sets[@]}"

# Set the number of threads
num_threads=4

# Check if the record file exists
if [ -f "record.txt" ]; then
  # Read the record file and determine the start point
  start_point=$(wc -l < "record.txt")
else
  start_point=1
fi

# Function to process a gene file
process_gene_file() {
  gene_file=$1
  gene_name=$(echo $gene_file | sed -E 's/(UCE-[0-9]+).*/\1/')
  gene_path="$geneset_loc$gene_file"
  store_path="./test/$gene_name/"

  # Set the path and filename for iqtree running result
  prefix_single="${store_path}Single_${gene_name}"
  prefix_mix="${store_path}Mix_${gene_name}"
  # Set the command for both single and mix model
  arg_single=("-s" "$gene_path" "-B" "1000" "--prefix" "$prefix_single")
  arg_mix=("-s" "$gene_path" "-m" "ESTMIXNUM" "-mrate" "E,I,G,I+G,R,I+R" "-opt_qmix_criteria" "1" "--prefix" "$prefix_mix")
  
  # Create the store_path directory if it doesn't exist
  mkdir -p "$store_path"

  # Command to run the single-class model in iqtree
  "$iqtree_loc" "${arg_single[@]}" > /dev/null
  # Command to run the mixture-class model in iqtree
  "$iqtree_loc" "${arg_mix[@]}" > /dev/null
  
  # Save the necessary information
  export outgroup species_tree_file gene_file gene_path gene_name prefix_single prefix_mix
  Rscript -e 'save(list = c(Sys.getenv("outgroup"), Sys.getenv("species_tree_file"), Sys.getenv("gene_file"), Sys.getenv("gene_path"), Sys.getenv("gene_name"), Sys.getenv("prefix_single"), Sys.getenv("prefix_mix")), file = paste0(Sys.getenv("store_path"), "path_info.RData"))'
  # Render the R Markdown file
  Rscript -e 'rmarkdown::render("tree_comparison_combined.Rmd", params = list(workingdict = Sys.getenv("store_path")), output_file = paste0(Sys.getenv("store_path"), Sys.getenv("gene_name"), "_summary.html"))'
  
  # Append the gene record to the record file
  echo "$gene_name" >> "record.txt"
}

# Export the function to be used by parallel
export -f process_gene_file

# Process gene files in parallel
echo "Processing gene files..."
parallel -j $num_threads process_gene_file ::: "${gene_sets[@]:$start_point-1}"
echo "Processing completed."
