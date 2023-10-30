#!/bin/bash

# Set the paths and variables
iqtree_loc="/home/tim/software/iqtree-2.2.2.7.modelmix-Linux/bin/iqtree2"
geneset_loc="/home/tim/mixmodel/datasets/different_method/Squamates-Final_Alns/MAFFT-Auto-untrimmed/"
outgroup_path="/home/tim/mixmodel/datasets/different_method/Squamates-Final_Alns/Squamates-Taxon-Outgroup.txt"
Species_tree_path="/home/tim/mixmodel/datasets/different_method/phylogenic_tree/Squamates/Rooted_RAxML_bipartitions.MAFFTAuto_Untrimmed.tre"

# Read the gene sets
gene_sets=( $(ls $geneset_loc | sort) )

echo "Current gene file location: $geneset_loc"
echo "Number of genes: ${#gene_sets[@]}"

# Set the number of threads
num_threads=20

# Set the execution count (the number of files to process)
execution_count=4430

# Set restart flag (True or False)
restart=True
skip_down=False

# Export the variables to be used by parallel
export iqtree_loc geneset_loc outgroup_path Species_tree_path restart skip_down

# If restart is True, delete the record.txt file if it exists
if [[ $restart = True && -f "record.txt" ]]; then
  rm "record.txt"
fi

# If the record file is not empty, read it and remove those genes from gene_sets
if [ -s "record.txt" ]; then
  mapfile -t processed_genes < record.txt
  for processed_gene in "${processed_genes[@]}"; do
    # Find the index of the processed gene in the gene_sets array
    for i in "${!gene_sets[@]}"; do
      if [[ ${gene_sets[i]} = $processed_gene ]]; then
        # Unset the element at the found index
        unset 'gene_sets[i]'
      fi
    done
  done
  # Re-index the array to remove gaps
  gene_sets=("${gene_sets[@]}")
fi

# Function to process a gene file
process_gene_file() {
  start_time=$(date +%s)
  gene_file=$1
  gene_name=$(echo $gene_file | sed -E 's/(UCE-[0-9]+).*/\1/')
  gene_path="$geneset_loc$gene_file"
  store_path="./result/$gene_name/"
  
  # Create the store_path directory if it doesn't exist
  mkdir -p "$store_path"

  # Check if HTML file exists and skip_down flag is False, then skip this gene file
  if [[ -f "${store_path}${gene_name}_summary.html" && $skip_down = True ]]; then
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
  printf "outgroup_path=%s\nSpecies_tree_path=%s\ngene_file=%s\ngene_path=%s\ngene_name=%s\nprefix_single=%s\nprefix_mix=%s\n" "$outgroup_path" "$Species_tree_path" "$gene_file" "$gene_path" "$gene_name" "$prefix_single" "$prefix_mix" > "${store_path}path_info.txt"
  
  # Render the R Markdown file
  Rscript -e "rmarkdown::render('tree_comparison_combined.Rmd', params = list(workingdict = '${store_path}'), output_file = paste0('${store_path}', '${gene_name}', '_summary.html'))" > "${store_path}${gene_name}_rmd.log" 2>&1
  time_diff=$((($(date +%s) - $start_time)/60))
  
  # Check if HTML file was created
  if [[ -f "${store_path}${gene_name}_summary.html" ]]; then
    # Append the gene record to the record file, use full filename instead of just gene name
    echo "$gene_file" >> "record.txt"
    echo "Completed processing for $gene_name in $time_diff minutes. Total genes processed so far: $(wc -l < 'record.txt')."
  else
    # Print the error log
    cat "${store_path}${gene_name}_rmd.log"
    echo "Failed processing for $gene_name."
  fi
}

# Export the function to be used by parallel
export -f process_gene_file

# Process gene files in parallel, only process up to execution_count files.
echo "Processing gene files..."
parallel --env process_gene_file -j $num_threads process_gene_file ::: "${gene_sets[@]:0:$execution_count}" 
echo "Processing completed."
