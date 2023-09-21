
# "geneset_loc" is the folder where you save a set of multi-alignments
iqtree_loc <- "/home/timmu/miniconda3/envs/phylo/bin/iqtree2"
geneset_loc <- "../datasets/different_method/Squamates-Final_Alns/MAFFT-Auto-gt/"
# The species tree file should be assigned manually
outgroup <- readLines("../datasets/different_method/Squamates-Final_Alns/Squamates-Taxon-Outgroup.txt")
Species_tree <- read.tree("/home/timmu/project/mixmodel/datasets/different_method/phylogenic_tree/Squamates/Rooted_MAFFT_Auto_GapThreshold_Species_Tree_Scored.tre")

gene_sets <- system2(command = "ls", args = geneset_loc, stdout = TRUE)
gene_sets <- sort(gene_sets)
cat("Current gene file location:", geneset_loc, "\n")
cat("Number of gene:", length(gene_sets), "\n")

# Add a parameter restart, if TRUE, then start from the beginning
restart <- FALSE

# Check if there is a save file
start_point <- 1
if (file.exists("checkpoint.RData") & !restart) {
  load("checkpoint.RData")
}

end_point <- 30

for (i in start_point:min(length(gene_sets), end_point)){
  candidate <- i
  id <- as.character(candidate)
  gene_file <- gene_sets[candidate]
  gene_path <- paste(geneset_loc, gene_file, sep = "")
  # Modify this to get gene_name
  gene_name <- gsub("(UCE-\\d+).*", "\\1", gene_file)
  cat("Current file name:", gene_sets[candidate], "\n")
  store_path <- paste("./test/", gene_name, "/", sep = "")
    
  # Set the path and filename for iqtree running result
  prefix_single <- paste(store_path, "Single_",gene_name, sep = "")
  prefix_mix <- paste(store_path, "Mix_", gene_name, sep = "")
  # Set the command for both single and mix model
  arg_single <- c("-s", gene_path, "-B", "1000", "--prefix", prefix_single)
  arg_mix <- c("-s", gene_path, "-m", "ESTMIXNUM", "-mrate", "E,I,G,I+G,R,I+R", "-opt_qmix_criteria", "1", "--prefix", prefix_mix)
  dir.create(store_path, showWarnings = FALSE)
  
  # Command of run one class model iq-tree
  system2(command = iqtree_loc, args = arg_single, stdout = FALSE)
  # Command of run mixture class model iq-tree
  system2(command = iqtree_loc, args = arg_mix, stdout = FALSE)
  
  # Save the value of candidate at the end of each loop
  save(outgroup, Species_tree, gene_file, gene_path, gene_name, prefix_single, prefix_mix, file = paste(store_path, "path_info.RData", sep = ""))
  if (file.exists(store_path)) {
    if (restart) {
      rmarkdown::render("tree_comparison_combined.Rmd", params = list(workingdict = store_path), output_file = paste(store_path, gene_name,"_summary.html", sep = ""))    } else {
    }
  } else {
    rmarkdown::render("tree_comparison_combined.Rmd", params = list(workingdict = store_path), output_file = paste(store_path, gene_name,"_summary.html", sep = ""))
  }

  save(i, file = "checkpoint.RData")
}


