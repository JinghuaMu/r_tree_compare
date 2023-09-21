#!/bin/bash

# Set the paths and variables
iqtree_loc="/home/timmu/miniconda3/envs/phylo/bin/iqtree2"
geneset_loc="../datasets/different_method/Squamates-Final_Alns/MAFFT-Auto-gt/"
outgroup_file="../datasets/different_method/Squamates-Final_Alns/Squamates-Taxon-Outgroup.txt"
Species_tree_path="/home/timmu/project/mixmodel/datasets/different_method/phylogenic_tree/Squamates/Rooted_MAFFT_Auto_GapThreshold_Species_Tree_Scored.tre"

processed_genes=0
# Export the variables to be used by parallel
export iqtree_loc geneset_loc outgroup_file Species_tree_path processed_genes

# Read the gene sets
gene_sets=( $(ls $geneset_loc | sort) )

echo "Current gene file location: $geneset_loc"
echo "Number of genes: ${#gene_sets[@]}"

# Set the number of threads
num_threads=10

# Check if the record file exists
if [ -f "record.txt" ]; then
  # Read the record file and determine the start point
  start_point=$(wc -l < "record.txt")
else
  start_point=1
fi

# Set the end point
end_point=31

# Set restart flag (True or False)
restart=False

# Function to process a gene file
process_gene_file() {
  gene_file=$1
  gene_name=$(echo $gene_file | sed -E 's/(UCE-[0-9]+).*/\1/')
  gene_path="$geneset_loc$gene_file"
  store_path="./test/$gene_name/"

  # Create the store_path directory if it doesn't exist
  mkdir -p "$store_path"

  # Check if HTML file exists and restart flag is False, then skip this gene file
  if [[ -f "${store_path}${gene_name}_summary.html" && $restart = False ]]; then
    echo "HTML file for $gene_name already exists, skipping..."
    return
  fi

  # Set the path and filename for iqtree running result
  prefix_single="${store_path}Single_${gene_name}"
  prefix_mix="${store_path}Mix_${gene_name}"
  
  # Set the command for both single and mix model
  arg_single=("-s" "$gene_path" "-B" "1000" "--prefix" "$prefix_single")
  arg_mix=("-s" "$gene_path" "-m" "ESTMIXNUM" "-mrate" "E,I,G,I+G,R,I+R" "-opt_qmix_criteria" "1" "--prefix" "$prefix_mix")

  # Command to run the single-class model in iqtree
  "$iqtree_loc" "${arg_single[@]}" &> /dev/null
  
  # Command to run the mixture-class model in iqtree
  "$iqtree_loc" "${arg_mix[@]}" &> /dev/null
  
  # Save the necessary information as a txt file instead of RData file
  printf "outgroup=%s\nSpecies_tree_path=%s\ngene_file=%s\ngene_path=%s\ngene_name=%s\nprefix_single=%s\nprefix_mix=%s\n" "$outgroup_file" "$Species_tree_path" "$gene_file" "$gene_path" "$gene_name" "$prefix_single" "$prefix_mix" > "${store_path}path_info.txt"
  
  # Render the R Markdown file
  Rscript -e "rmarkdown::render('tree_comparison_combined.Rmd', params = list(workingdict = '${store_path}'), output_file = paste0('${store_path}', '${gene_name}', '_summary.html'))" > "${store_path}${gene_name}_rmd.log" 2>&1
  
  # Append the gene record to the record file
  echo "$gene_name" >> "record.txt"

  # Increment processed gene counter and print a message showing completion of this gene and total completed so far
  ((processed_genes++))
  echo "Completed processing for $gene_name. Total genes processed so far: $processed_genes."
}

# Export the function to be used by parallel
export -f process_gene_file

# Process gene files in parallel, only process up to end_point files.
echo "Processing gene files..."
parallel --env process_gene_file -j $num_threads process_gene_file ::: "${gene_sets[@]:$start_point-1:$end_point-1}"
echo "Processing completed."
